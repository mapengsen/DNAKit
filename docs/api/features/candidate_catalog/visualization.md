# Visualization

This category contains **6** documented candidate capabilities; 6 have a locally callable formal entry (0 partial due to binary-format dependencies; 0 bounded external backends; 6 dependency-adapter only). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) export_figure

- **Purpose:** Document a reusable export figure capability for Visualization workflows. Candidate CAP-0822 comes from 45 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the export figure operation, and return structured results with provenance and diagnostics. The optional `plot_kind` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/visualization/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.visualization.export_figure(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `export_figure` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.visualization import export_figure

result = export_figure(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0822",
  "operation": "export_figure",
  "result_type": "dict",
  "result_preview": "{'module': 'visualization', 'operation': 'export_figure', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 2) generate_analysis_summary

- **Purpose:** Document a reusable generate analysis summary capability for Visualization workflows. Candidate CAP-0823 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the generate analysis summary operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/visualization/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.visualization.generate_analysis_summary(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `generate_analysis_summary` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.visualization import generate_analysis_summary

result = generate_analysis_summary(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0823",
  "operation": "generate_analysis_summary",
  "result_type": "dict",
  "result_preview": "{'module': 'visualization', 'operation': 'generate_analysis_summary', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7c...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 3) generate_clinical_pdf_report

- **Purpose:** Document a reusable generate clinical pdf report capability for Visualization workflows. Candidate CAP-0824 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the generate clinical pdf report operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/visualization/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.visualization.generate_clinical_pdf_report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `generate_clinical_pdf_report` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.visualization import generate_clinical_pdf_report

result = generate_clinical_pdf_report(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0824",
  "operation": "generate_clinical_pdf_report",
  "result_type": "dict",
  "result_preview": "{'module': 'visualization', 'operation': 'generate_clinical_pdf_report', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49d...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 4) generate_html_genome_report

- **Purpose:** Document a reusable generate html genome report capability for Visualization workflows. Candidate CAP-0825 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the generate html genome report operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/visualization/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.visualization.generate_html_genome_report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `generate_html_genome_report` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.visualization import generate_html_genome_report

result = generate_html_genome_report(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0825",
  "operation": "generate_html_genome_report",
  "result_type": "dict",
  "result_preview": "{'module': 'visualization', 'operation': 'generate_html_genome_report', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 5) generate_qc_html_report

- **Purpose:** Document a reusable generate qc html report capability for Visualization workflows. Candidate CAP-0826 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the generate qc html report operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/visualization/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.visualization.generate_qc_html_report(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `generate_qc_html_report` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.visualization import generate_qc_html_report

result = generate_qc_html_report(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0826",
  "operation": "generate_qc_html_report",
  "result_type": "dict",
  "result_preview": "{'module': 'visualization', 'operation': 'generate_qc_html_report', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 6) save_figure

- **Purpose:** Document a reusable save figure capability for Visualization workflows. Candidate CAP-0828 comes from 180 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the save figure operation, and return structured results with provenance and diagnostics. The optional `plot_kind` parameter selects one supported variant (the supported values recorded in the source catalog). The formal package reuses the existing Matplotlib figure persistence API; plot-kind and output-format variants remain options of the existing callable. Review basis: `generic_operation_name`. Existing-name evidence: `dnakit.visualization.save_figure|save_figure`.

- **API:** `dnakit.visualization.save_figure(figure[required], target[required], overwrite[optional], tight_layout[optional], close[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `figure`, `target`; optional arguments: `overwrite`, `tight_layout`, `close`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The existing formal `save_figure` function returns the staged output `Path`. It accepts a Matplotlib figure and chooses the output format from the target suffix.

- **Example code:**

```python
from pathlib import Path

import matplotlib.pyplot as plt

from dnakit.visualization import save_figure

figure, axis = plt.subplots()
axis.plot([0, 1], [0, 1])
result = save_figure(figure, Path("results/save_figure.svg"), overwrite=True)
print(result)
```

- **Example output:**

```json
{
  "status": "formal_existing",
  "capability_id": "CAP-0828",
  "operation": "save_figure",
  "result_type": "PosixPath",
  "result_preview": "results/save_figure.svg",
  "provenance": {
    "schema": "formal_dnakit",
    "implementation": "existing_visualization_api"
  }
}
```
