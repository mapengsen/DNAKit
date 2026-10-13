# Sequencing quality control

This category contains **50** documented candidate capabilities; 50 have a locally callable formal entry (0 partial due to binary-format dependencies; 4 bounded external backends; 46 dependency-adapter only). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) assign_sample_barcodes

- **Purpose:** Document a reusable assign sample barcodes capability for Sequencing quality control workflows. Candidate CAP-0637 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the assign sample barcodes operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.assign_sample_barcodes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `assign_sample_barcodes` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import assign_sample_barcodes

result = assign_sample_barcodes(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0637",
  "operation": "assign_sample_barcodes",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'assign_sample_barcodes', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 2) basecall_nanopore_reads

- **Purpose:** Document a reusable basecall nanopore reads capability for Sequencing quality control workflows. Candidate CAP-0638 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the basecall nanopore reads operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.basecall_nanopore_reads(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `basecall_nanopore_reads` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import basecall_nanopore_reads

result = basecall_nanopore_reads(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0638",
  "operation": "basecall_nanopore_reads",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'basecall_nanopore_reads', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 3) calculate_adapter_content

- **Purpose:** Document a reusable calculate adapter content capability for Sequencing quality control workflows. Candidate CAP-0639 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate adapter content operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_adapter_content(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_adapter_content` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_adapter_content

result = calculate_adapter_content(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0639",
  "operation": "calculate_adapter_content",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_adapter_content', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7c...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 4) calculate_adapter_content_tenx

- **Purpose:** Document a reusable calculate adapter content tenx capability for Sequencing quality control workflows. Candidate CAP-0640 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate adapter content tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_adapter_content_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_adapter_content_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_adapter_content_tenx

result = calculate_adapter_content_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0640",
  "operation": "calculate_adapter_content_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_adapter_content_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a4...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 5) calculate_barcode_quality

- **Purpose:** Document a reusable calculate barcode quality capability for Sequencing quality control workflows. Candidate CAP-0641 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate barcode quality operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_barcode_quality(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_barcode_quality` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_barcode_quality

result = calculate_barcode_quality(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0641",
  "operation": "calculate_barcode_quality",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_barcode_quality', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7c...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 6) calculate_barcode_quality_tenx

- **Purpose:** Document a reusable calculate barcode quality tenx capability for Sequencing quality control workflows. Candidate CAP-0642 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate barcode quality tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_barcode_quality_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_barcode_quality_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_barcode_quality_tenx

result = calculate_barcode_quality_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0642",
  "operation": "calculate_barcode_quality_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_barcode_quality_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a4...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 7) calculate_basecalling_qc

- **Purpose:** Document a reusable calculate basecalling qc capability for Sequencing quality control workflows. Candidate CAP-0643 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate basecalling qc operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_basecalling_qc(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_basecalling_qc` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_basecalling_qc

result = calculate_basecalling_qc(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0643",
  "operation": "calculate_basecalling_qc",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_basecalling_qc', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 8) calculate_basecalling_qc_tenx

- **Purpose:** Document a reusable calculate basecalling qc tenx capability for Sequencing quality control workflows. Candidate CAP-0644 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate basecalling qc tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_basecalling_qc_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_basecalling_qc_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_basecalling_qc_tenx

result = calculate_basecalling_qc_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0644",
  "operation": "calculate_basecalling_qc_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_basecalling_qc_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 9) calculate_contamination

- **Purpose:** Document a reusable calculate contamination capability for Sequencing quality control workflows. Candidate CAP-0645 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate contamination operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_contamination(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_contamination` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_contamination

result = calculate_contamination(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0645",
  "operation": "calculate_contamination",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_contamination', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 10) calculate_contamination_tenx

- **Purpose:** Document a reusable calculate contamination tenx capability for Sequencing quality control workflows. Candidate CAP-0646 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate contamination tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_contamination_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_contamination_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_contamination_tenx

result = calculate_contamination_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0646",
  "operation": "calculate_contamination_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_contamination_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49d...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 11) calculate_duplicate_rate

- **Purpose:** Document a reusable calculate duplicate rate capability for Sequencing quality control workflows. Candidate CAP-0647 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate duplicate rate operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_duplicate_rate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_duplicate_rate` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_duplicate_rate

result = calculate_duplicate_rate(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0647",
  "operation": "calculate_duplicate_rate",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_duplicate_rate', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 12) calculate_duplication

- **Purpose:** Document a reusable calculate duplication capability for Sequencing quality control workflows. Candidate CAP-0648 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate duplication operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_duplication(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_duplication` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_duplication

result = calculate_duplication(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0648",
  "operation": "calculate_duplication",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_duplication', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd60...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 13) calculate_duplication_tenx

- **Purpose:** Document a reusable calculate duplication tenx capability for Sequencing quality control workflows. Candidate CAP-0649 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate duplication tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_duplication_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_duplication_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_duplication_tenx

result = calculate_duplication_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0649",
  "operation": "calculate_duplication_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_duplication_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 14) calculate_error_profile

- **Purpose:** Document a reusable calculate error profile capability for Sequencing quality control workflows. Candidate CAP-0650 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate error profile operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_error_profile(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_error_profile` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_error_profile

result = calculate_error_profile(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0650",
  "operation": "calculate_error_profile",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_error_profile', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 15) calculate_error_profile_tenx

- **Purpose:** Document a reusable calculate error profile tenx capability for Sequencing quality control workflows. Candidate CAP-0651 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate error profile tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_error_profile_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_error_profile_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_error_profile_tenx

result = calculate_error_profile_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0651",
  "operation": "calculate_error_profile_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_error_profile_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49d...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 16) calculate_gc_bias

- **Purpose:** Document a reusable calculate gc bias capability for Sequencing quality control workflows. Candidate CAP-0652 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate gc bias operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_gc_bias(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_gc_bias` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_gc_bias

result = calculate_gc_bias(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0652",
  "operation": "calculate_gc_bias",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_gc_bias', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f8...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 17) calculate_gc_bias_tenx

- **Purpose:** Document a reusable calculate gc bias tenx capability for Sequencing quality control workflows. Candidate CAP-0653 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate gc bias tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_gc_bias_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_gc_bias_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_gc_bias_tenx

result = calculate_gc_bias_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0653",
  "operation": "calculate_gc_bias_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_gc_bias_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 18) calculate_index_hopping

- **Purpose:** Document a reusable calculate index hopping capability for Sequencing quality control workflows. Candidate CAP-0654 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate index hopping operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_index_hopping(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_index_hopping` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_index_hopping

result = calculate_index_hopping(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0654",
  "operation": "calculate_index_hopping",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_index_hopping', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 19) calculate_index_hopping_tenx

- **Purpose:** Document a reusable calculate index hopping tenx capability for Sequencing quality control workflows. Candidate CAP-0655 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate index hopping tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_index_hopping_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_index_hopping_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_index_hopping_tenx

result = calculate_index_hopping_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0655",
  "operation": "calculate_index_hopping_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_index_hopping_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49d...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 20) calculate_insert_size

- **Purpose:** Document a reusable calculate insert size capability for Sequencing quality control workflows. Candidate CAP-0656 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate insert size operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_insert_size(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_insert_size` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_insert_size

result = calculate_insert_size(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0656",
  "operation": "calculate_insert_size",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_insert_size', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd60...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 21) calculate_insert_size_distribution

- **Purpose:** Document a reusable calculate insert size distribution capability for Sequencing quality control workflows. Candidate CAP-0657 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate insert size distribution operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_insert_size_distribution(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_insert_size_distribution` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_insert_size_distribution

result = calculate_insert_size_distribution(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0657",
  "operation": "calculate_insert_size_distribution",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_insert_size_distribution', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 22) calculate_insert_size_tenx

- **Purpose:** Document a reusable calculate insert size tenx capability for Sequencing quality control workflows. Candidate CAP-0658 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate insert size tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_insert_size_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_insert_size_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_insert_size_tenx

result = calculate_insert_size_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0658",
  "operation": "calculate_insert_size_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_insert_size_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 23) calculate_phasing_qc

- **Purpose:** Document a reusable calculate phasing qc capability for Sequencing quality control workflows. Candidate CAP-0659 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate phasing qc operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_phasing_qc(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_phasing_qc` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_phasing_qc

result = calculate_phasing_qc(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0659",
  "operation": "calculate_phasing_qc",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_phasing_qc', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd602...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 24) calculate_phasing_qc_tenx

- **Purpose:** Document a reusable calculate phasing qc tenx capability for Sequencing quality control workflows. Candidate CAP-0660 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate phasing qc tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_phasing_qc_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_phasing_qc_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_phasing_qc_tenx

result = calculate_phasing_qc_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0660",
  "operation": "calculate_phasing_qc_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_phasing_qc_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7c...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 25) calculate_qscore_summary

- **Purpose:** Document a reusable calculate qscore summary capability for Sequencing quality control workflows. Candidate CAP-0661 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate qscore summary operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_qscore_summary(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_qscore_summary` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_qscore_summary

result = calculate_qscore_summary(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0661",
  "operation": "calculate_qscore_summary",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_qscore_summary', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 26) calculate_quality

- **Purpose:** Document a reusable calculate quality capability for Sequencing quality control workflows. Candidate CAP-0662 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate quality operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. Existing-name evidence: `evaluate_quality`.

- **API:** `dnakit.sequencing_qc.calculate_quality(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_quality` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_quality

result = calculate_quality(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0662",
  "operation": "calculate_quality",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_quality', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f8...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 27) calculate_quality_tenx

- **Purpose:** Document a reusable calculate quality tenx capability for Sequencing quality control workflows. Candidate CAP-0663 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate quality tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_quality_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_quality_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_quality_tenx

result = calculate_quality_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0663",
  "operation": "calculate_quality_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_quality_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 28) calculate_read_length

- **Purpose:** Document a reusable calculate read length capability for Sequencing quality control workflows. Candidate CAP-0664 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate read length operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_read_length(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_read_length` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_read_length

result = calculate_read_length(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0664",
  "operation": "calculate_read_length",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_read_length', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd60...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 29) calculate_read_length_tenx

- **Purpose:** Document a reusable calculate read length tenx capability for Sequencing quality control workflows. Candidate CAP-0665 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate read length tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_read_length_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_read_length_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_read_length_tenx

result = calculate_read_length_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0665",
  "operation": "calculate_read_length_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_read_length_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 30) calculate_yield

- **Purpose:** Document a reusable calculate yield capability for Sequencing quality control workflows. Candidate CAP-0666 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate yield operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_yield(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_yield` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_yield

result = calculate_yield(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0666",
  "operation": "calculate_yield",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_yield', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 31) calculate_yield_tenx

- **Purpose:** Document a reusable calculate yield tenx capability for Sequencing quality control workflows. Candidate CAP-0667 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate yield tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.calculate_yield_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_yield_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import calculate_yield_tenx

result = calculate_yield_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0667",
  "operation": "calculate_yield_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'calculate_yield_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd602...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 32) collapse_umi_reads

- **Purpose:** Document a reusable collapse umi reads capability for Sequencing quality control workflows. Candidate CAP-0668 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the collapse umi reads operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.collapse_umi_reads(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `collapse_umi_reads` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import collapse_umi_reads

result = collapse_umi_reads(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0668",
  "operation": "collapse_umi_reads",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'collapse_umi_reads', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 33) correct_long_reads

- **Purpose:** Document a reusable correct long reads capability for Sequencing quality control workflows. Candidate CAP-0669 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the correct long reads operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.correct_long_reads(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `correct_long_reads` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import correct_long_reads

result = correct_long_reads(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0669",
  "operation": "correct_long_reads",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'correct_long_reads', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 34) correct_umi_errors

- **Purpose:** Document a reusable correct umi errors capability for Sequencing quality control workflows. Candidate CAP-0670 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the correct umi errors operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.correct_umi_errors(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `correct_umi_errors` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import correct_umi_errors

result = correct_umi_errors(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0670",
  "operation": "correct_umi_errors",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'correct_umi_errors', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 35) demultiplex_barcoded_reads

- **Purpose:** Document a reusable demultiplex barcoded reads capability for Sequencing quality control workflows. Candidate CAP-0671 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the demultiplex barcoded reads operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.demultiplex_barcoded_reads(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `demultiplex_barcoded_reads` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import demultiplex_barcoded_reads

result = demultiplex_barcoded_reads(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0671",
  "operation": "demultiplex_barcoded_reads",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'demultiplex_barcoded_reads', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 36) detect_cross_sample_contamination

- **Purpose:** Document a reusable detect cross sample contamination capability for Sequencing quality control workflows. Candidate CAP-0672 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect cross sample contamination operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.detect_cross_sample_contamination(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `detect_cross_sample_contamination` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import detect_cross_sample_contamination

result = detect_cross_sample_contamination(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0672",
  "operation": "detect_cross_sample_contamination",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'detect_cross_sample_contamination', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b4...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 37) detect_platform

- **Purpose:** Document a reusable detect platform capability for Sequencing quality control workflows. Candidate CAP-0673 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect platform operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.detect_platform(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `detect_platform` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import detect_platform

result = detect_platform(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0673",
  "operation": "detect_platform",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'detect_platform', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 38) detect_platform_tenx

- **Purpose:** Document a reusable detect platform tenx capability for Sequencing quality control workflows. Candidate CAP-0674 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect platform tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.detect_platform_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `detect_platform_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import detect_platform_tenx

result = detect_platform_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0674",
  "operation": "detect_platform_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'detect_platform_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd602...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 39) estimate_library_complexity

- **Purpose:** Document a reusable estimate library complexity capability for Sequencing quality control workflows. Candidate CAP-0675 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the estimate library complexity operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.estimate_library_complexity(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `estimate_library_complexity` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import estimate_library_complexity

result = estimate_library_complexity(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0675",
  "operation": "estimate_library_complexity",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'estimate_library_complexity', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 40) estimate_read_error_model

- **Purpose:** Document a reusable estimate read error model capability for Sequencing quality control workflows. Candidate CAP-0676 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the estimate read error model operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.estimate_read_error_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `estimate_read_error_model` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import estimate_read_error_model

result = estimate_read_error_model(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0676",
  "operation": "estimate_read_error_model",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'estimate_read_error_model', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7c...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 41) extract_umi_sequences

- **Purpose:** Document a reusable extract umi sequences capability for Sequencing quality control workflows. Candidate CAP-0677 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract umi sequences operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.extract_umi_sequences(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `extract_umi_sequences` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import extract_umi_sequences

result = extract_umi_sequences(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0677",
  "operation": "extract_umi_sequences",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'extract_umi_sequences', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd60...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 42) generate_qc_report

- **Purpose:** Document a reusable generate qc report capability for Sequencing quality control workflows. Candidate CAP-0678 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the generate qc report operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.generate_qc_report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `generate_qc_report` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import generate_qc_report

result = generate_qc_report(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0678",
  "operation": "generate_qc_report",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'generate_qc_report', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 43) generate_qc_report_tenx

- **Purpose:** Document a reusable generate qc report tenx capability for Sequencing quality control workflows. Candidate CAP-0679 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the generate qc report tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.generate_qc_report_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `generate_qc_report_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import generate_qc_report_tenx

result = generate_qc_report_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0679",
  "operation": "generate_qc_report_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'generate_qc_report_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 44) plot_per_base_quality

- **Purpose:** Document a reusable plot per base quality capability for Sequencing quality control workflows. Candidate CAP-0680 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the plot per base quality operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.plot_per_base_quality(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `plot_per_base_quality` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import plot_per_base_quality

result = plot_per_base_quality(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0680",
  "operation": "plot_per_base_quality",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'plot_per_base_quality', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd60...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 45) plot_qc_summary

- **Purpose:** Document a reusable plot qc summary capability for Sequencing quality control workflows. Candidate CAP-0681 comes from 13 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the plot qc summary operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.plot_qc_summary(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `plot_qc_summary` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import plot_qc_summary

result = plot_qc_summary(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0681",
  "operation": "plot_qc_summary",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'plot_qc_summary', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 46) plot_qc_summary_tenx

- **Purpose:** Document a reusable plot qc summary tenx capability for Sequencing quality control workflows. Candidate CAP-0682 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the plot qc summary tenx operation, and return structured results with provenance and diagnostics. The optional `platform` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/sequencing_qc/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.sequencing_qc.plot_qc_summary_tenx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `plot_qc_summary_tenx` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import plot_qc_summary_tenx

result = plot_qc_summary_tenx(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0682",
  "operation": "plot_qc_summary_tenx",
  "result_type": "dict",
  "result_preview": "{'module': 'sequencing_qc', 'operation': 'plot_qc_summary_tenx', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd602...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 47) run_cutadapt

- **Purpose:** Document a reusable run cutadapt capability for Sequencing quality control workflows. Candidate CAP-0683 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the local paths and typed options, invoke the caller-installed Cutadapt executable without a shell, enforce time/output limits, verify declared trimmed FASTQ files and an optional JSON report, and return structured provenance. The fixed option allowlist is a DNAKit boundary rather than a claim to expose the backend's complete CLI. The wrapper records the Martin (2011), DOI 10.14806/ej.17.1.200 citation and MIT license evidence.

- **API:** `dnakit.sequencing_qc.run_cutadapt(input, *, options=None, verbose=True)` (formal DNAKit namespace; the same bounded external-process effect is registered for Agent/MCP).

- **Input:** Required `input`: a FASTQ path or a `{read1, read2}` mapping. Required `options` values include the caller-installed `executable_path` and output location(s); `verbose` is optional. Only the documented allow-listed options are accepted, and output paths must not already exist.

- **Output:** The wrapper returns a JSON-serializable completion envelope with `status`, `backend`, the exact argv tokens, declared `outputs`, bounded backend output, and `provenance`; the backend produces declared trimmed FASTQ files and an optional JSON report. Missing executables, non-zero exits, missing artifacts, invalid JSON reports and output-limit violations are reported as structured errors.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import run_cutadapt

result = run_cutadapt(
    input={"read1": Path("reads_R1.fastq"), "read2": Path("reads_R2.fastq")},
    options={
        "executable_path": "/usr/local/bin/cutadapt",
        "output1": Path("results/trimmed_R1.fastq"),
        "output2": Path("results/trimmed_R2.fastq"),
        "adapters": ["AGATCGGAAG"],
        "json_report": Path("results/cutadapt.json"),
    },
    verbose=False,
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "backend": "cutadapt",
  "executable": "/usr/local/bin/backend",
  "arguments": [
    "<validated argv tokens>"
  ],
  "outputs": [
    "results/report.html"
  ],
  "provenance": {
    "implementation": "bounded_external_adapter",
    "citation": "Martin (2011), DOI 10.14806/ej.17.1.200",
    "license": "MIT"
  }
}
```

## 48) run_fastp

- **Purpose:** Document a reusable run fastp capability for Sequencing quality control workflows. Candidate CAP-0684 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the local paths and typed options, invoke the caller-installed fastp executable without a shell, enforce time/output limits, verify declared filtered FASTQ files plus JSON and HTML reports, and return structured provenance. The fixed option allowlist is a DNAKit boundary rather than a claim to expose the backend's complete CLI. The wrapper records the Chen (2025), DOI 10.1002/imt2.70078 citation and MIT license evidence.

- **API:** `dnakit.sequencing_qc.run_fastp(input, *, options=None, verbose=True)` (formal DNAKit namespace; the same bounded external-process effect is registered for Agent/MCP).

- **Input:** Required `input`: a FASTQ path or a `{read1, read2}` mapping. Required `options` values include the caller-installed `executable_path` and output location(s); `verbose` is optional. Only the documented allow-listed options are accepted, and output paths must not already exist.

- **Output:** The wrapper returns a JSON-serializable completion envelope with `status`, `backend`, the exact argv tokens, declared `outputs`, bounded backend output, and `provenance`; the backend produces declared filtered FASTQ files plus JSON and HTML reports. Missing executables, non-zero exits, missing artifacts, invalid JSON reports and output-limit violations are reported as structured errors.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import run_fastp

result = run_fastp(
    input=Path("reads.fastq"),
    options={
        "executable_path": "/usr/local/bin/fastp",
        "output": Path("results/filtered.fastq"),
        "json_report": Path("results/fastp.json"),
        "html_report": Path("results/fastp.html"),
    },
    verbose=False,
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "backend": "fastp",
  "executable": "/usr/local/bin/backend",
  "arguments": [
    "<validated argv tokens>"
  ],
  "outputs": [
    "results/report.html"
  ],
  "provenance": {
    "implementation": "bounded_external_adapter",
    "citation": "Chen (2025), DOI 10.1002/imt2.70078",
    "license": "MIT"
  }
}
```

## 49) run_fastqc

- **Purpose:** Document a reusable run fastqc capability for Sequencing quality control workflows. Candidate CAP-0685 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the local paths and typed options, invoke the caller-installed FastQC 0.12.x executable without a shell, enforce time/output limits, verify FastQC HTML and ZIP reports in an existing output directory, and return structured provenance. The fixed option allowlist is a DNAKit boundary rather than a claim to expose the backend's complete CLI. The wrapper records the FastQC project documentation; no DOI supplied by the survey citation and GPL-3.0 license evidence.

- **API:** `dnakit.sequencing_qc.run_fastqc(input, *, options=None, verbose=True)` (formal DNAKit namespace; the same bounded external-process effect is registered for Agent/MCP).

- **Input:** Required `input`: one to eight local FASTQ/BAM/SAM-style read paths. Required `options` values include the caller-installed `executable_path` and output location(s); `verbose` is optional. Only the documented allow-listed options are accepted, and output paths must not already exist.

- **Output:** The wrapper returns a JSON-serializable completion envelope with `status`, `backend`, the exact argv tokens, declared `outputs`, bounded backend output, and `provenance`; the backend produces FastQC HTML and ZIP reports in an existing output directory. Missing executables, non-zero exits, missing artifacts, invalid JSON reports and output-limit violations are reported as structured errors.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import run_fastqc

result = run_fastqc(
    input=Path("reads.fastq.gz"),
    options={
        "executable_path": "/opt/FastQC/fastqc",
        "output_directory": Path("results/fastqc"),
        "kmer": 7,
    },
    verbose=False,
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "backend": "fastqc",
  "executable": "/usr/local/bin/backend",
  "arguments": [
    "<validated argv tokens>"
  ],
  "outputs": [
    "results/report.html"
  ],
  "provenance": {
    "implementation": "bounded_external_adapter",
    "citation": "FastQC project documentation; no DOI supplied by the survey",
    "license": "GPL-3.0"
  }
}
```

## 50) run_multiqc

- **Purpose:** Document a reusable run multiqc capability for Sequencing quality control workflows. Candidate CAP-0686 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the local paths and typed options, invoke the caller-installed MultiQC executable without a shell, enforce time/output limits, verify a declared MultiQC HTML report in an empty output directory, and return structured provenance. The fixed option allowlist is a DNAKit boundary rather than a claim to expose the backend's complete CLI. The wrapper records the Ewels et al. (2016), DOI 10.1093/bioinformatics/btw354 citation and GPL-3.0-or-later license evidence.

- **API:** `dnakit.sequencing_qc.run_multiqc(input, *, options=None, verbose=True)` (formal DNAKit namespace; the same bounded external-process effect is registered for Agent/MCP).

- **Input:** Required `input`: one or more local files/directories containing QC inputs. Required `options` values include the caller-installed `executable_path` and output location(s); `verbose` is optional. Only the documented allow-listed options are accepted, and output paths must not already exist.

- **Output:** The wrapper returns a JSON-serializable completion envelope with `status`, `backend`, the exact argv tokens, declared `outputs`, bounded backend output, and `provenance`; the backend produces a declared MultiQC HTML report in an empty output directory. Missing executables, non-zero exits, missing artifacts, invalid JSON reports and output-limit violations are reported as structured errors.

- **Example code:**

```python
from pathlib import Path

from dnakit.sequencing_qc import run_multiqc

result = run_multiqc(
    input=[Path("results/fastqc")],
    options={
        "executable_path": "/opt/multiqc/bin/multiqc",
        "output_directory": Path("results/multiqc"),
        "filename": "multiqc_report.html",
    },
    verbose=False,
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "backend": "multiqc",
  "executable": "/usr/local/bin/backend",
  "arguments": [
    "<validated argv tokens>"
  ],
  "outputs": [
    "results/report.html"
  ],
  "provenance": {
    "implementation": "bounded_external_adapter",
    "citation": "Ewels et al. (2016), DOI 10.1093/bioinformatics/btw354",
    "license": "GPL-3.0-or-later"
  }
}
```
