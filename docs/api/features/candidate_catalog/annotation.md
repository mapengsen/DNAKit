# Annotation

This category contains **18** documented candidate capabilities; 18 have a locally callable formal entry (0 partial due to binary-format dependencies; 0 bounded external backends; 14 dependency-adapter only). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) annotate

- **Purpose:** Document a reusable annotate capability for Annotation workflows. Candidate CAP-0035 comes from 29 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the annotate operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (annotation_evidence, cds, centromere, enhancer, exon, gene, gene_family, insulator, intron, lncrna, mirna, orf, and 17 more). The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `annotate` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import annotate

result = annotate(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0035",
  "operation": "annotate",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'annotate', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed4389...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 2) assign_gene_ontology

- **Purpose:** Document a reusable assign gene ontology capability for Annotation workflows. Candidate CAP-0036 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the assign gene ontology operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.assign_gene_ontology(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `assign_gene_ontology` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import assign_gene_ontology

result = assign_gene_ontology(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0036",
  "operation": "assign_gene_ontology",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'assign_gene_ontology', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f8...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 3) assign_pathways

- **Purpose:** Document a reusable assign pathways capability for Annotation workflows. Candidate CAP-0037 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the assign pathways operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.assign_pathways(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `assign_pathways` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import assign_pathways

result = assign_pathways(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0037",
  "operation": "assign_pathways",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'assign_pathways', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 4) assign_protein_domains

- **Purpose:** Document a reusable assign protein domains capability for Annotation workflows. Candidate CAP-0038 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the assign protein domains operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.assign_protein_domains(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `assign_protein_domains` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import assign_protein_domains

result = assign_protein_domains(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0038",
  "operation": "assign_protein_domains",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'assign_protein_domains', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 5) build_gene_model

- **Purpose:** Document a reusable build gene model capability for Annotation workflows. Candidate CAP-0039 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the build gene model operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.build_gene_model(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `build_gene_model` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import build_gene_model

result = build_gene_model(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0039",
  "operation": "build_gene_model",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'build_gene_model', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 6) check_annotation_consistency

- **Purpose:** Document a reusable check annotation consistency capability for Annotation workflows. Candidate CAP-0040 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the check annotation consistency operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.check_annotation_consistency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `check_annotation_consistency` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import check_annotation_consistency

result = check_annotation_consistency(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0040",
  "operation": "check_annotation_consistency",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'check_annotation_consistency', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7c...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 7) compare

- **Purpose:** Document a reusable compare capability for Annotation workflows. Candidate CAP-0041 comes from 29 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compare operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (cds, centromere, enhancer, exon, gene, gene_family, gene_models, insulator, intron, lncrna, mirna, orf, and 17 more). The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `generic_operation_name`. Existing-name evidence: `compare|dnakit.similarity.compare`.

- **API:** `dnakit.annotation.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compare` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import compare

result = compare(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0041",
  "operation": "compare",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'compare', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 8) convert_annotation_coordinates

- **Purpose:** Document a reusable convert annotation coordinates capability for Annotation workflows. Candidate CAP-0042 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert annotation coordinates operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/annotation/local.py` and is exported through the formal DNAKit namespace. Coordinate conventions, contig lengths, and merge rules are explicit options; no annotation database or reference genome is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.convert_annotation_coordinates(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `convert_annotation_coordinates` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit coordinate provenance. It performs deterministic interval arithmetic only; no annotation database or reference sequence is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import convert_annotation_coordinates

result = convert_annotation_coordinates(
    input={"records": [{"chrom": "chr1", "start": 2, "end": 5, "id": "g1"}]},
    options={"coordinate_system": "1-based-closed", "target_coordinate_system": "0-based-half-open"}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0042",
  "operation": "convert_annotation_coordinates",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'convert_annotation_coordinates', 'result': {'records': [{'chrom': 'chr1', 'start': 1, 'end': 5, 'id': 'g1', 'coordinate_system': '0-based-half-open'}], 'record_count': 1, 'source_coordinate_system': '1-based-closed', 'target_coordinate_system': '0-based-half-...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 9) export

- **Purpose:** Document a reusable export capability for Annotation workflows. Candidate CAP-0043 comes from 28 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the export operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (cds, centromere, enhancer, exon, gene, gene_family, insulator, intron, lncrna, mirna, orf, ortholog, and 16 more). The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `export` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import export

result = export(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0043",
  "operation": "export",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'export', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed438948...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 10) extract

- **Purpose:** Document a reusable extract capability for Annotation workflows. Candidate CAP-0044 comes from 33 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (cds, cds_sequence, centromere, enhancer, exon, flanking_sequence, gene, gene_family, gene_sequence, insulator, intron, lncrna, and 21 more). The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.extract(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `extract` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import extract

result = extract(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0044",
  "operation": "extract",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'extract', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 11) find_homologous_genes

- **Purpose:** Document a reusable find homologous genes capability for Annotation workflows. Candidate CAP-0045 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the find homologous genes operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.find_homologous_genes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `find_homologous_genes` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import find_homologous_genes

result = find_homologous_genes(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0045",
  "operation": "find_homologous_genes",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'find_homologous_genes', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 12) index

- **Purpose:** Document a reusable index capability for Annotation workflows. Candidate CAP-0046 comes from 28 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the index operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (cds, centromere, enhancer, exon, gene, gene_family, insulator, intron, lncrna, mirna, orf, ortholog, and 16 more). The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.index(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `index` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import index

result = index(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0046",
  "operation": "index",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'index', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed4389481...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 13) make_intergenic_regions

- **Purpose:** Document a reusable make intergenic regions capability for Annotation workflows. Candidate CAP-0047 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the make intergenic regions operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/annotation/local.py` and is exported through the formal DNAKit namespace. Coordinate conventions, contig lengths, and merge rules are explicit options; no annotation database or reference genome is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.make_intergenic_regions(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `make_intergenic_regions` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit coordinate provenance. It performs deterministic interval arithmetic only; no annotation database or reference sequence is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import make_intergenic_regions

result = make_intergenic_regions(
    input={"genes": [{"chrom": "chr1", "start": 2, "end": 4}], "contig_lengths": {"chr1": 10}},
    options={"coordinate_system": "1-based-closed"}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0047",
  "operation": "make_intergenic_regions",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'make_intergenic_regions', 'result': {'regions': [{'id': 'intergenic_chr1_1', 'chrom': 'chr1', 'start': 1, 'end': 1, 'feature_type': 'intergenic', 'coordinate_system': '1-based-closed'}, {'id': 'intergenic_chr1_2', 'chrom': 'chr1', 'start': 5, 'end': 10, 'feat...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 14) make_promoter_regions

- **Purpose:** Document a reusable make promoter regions capability for Annotation workflows. Candidate CAP-0048 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the make promoter regions operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/annotation/local.py` and is exported through the formal DNAKit namespace. Coordinate conventions, contig lengths, and merge rules are explicit options; no annotation database or reference genome is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.make_promoter_regions(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `make_promoter_regions` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit coordinate provenance. It performs deterministic interval arithmetic only; no annotation database or reference sequence is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import make_promoter_regions

result = make_promoter_regions(
    input={"genes": [{"id": "g1", "chrom": "chr1", "start": 3, "end": 6, "strand": "+"}], "contig_lengths": {"chr1": 10}},
    options={"coordinate_system": "1-based-closed", "upstream": 2, "downstream": 1}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0048",
  "operation": "make_promoter_regions",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'make_promoter_regions', 'result': {'promoters': [{'id': 'promoter_g1', 'gene_id': 'g1', 'chrom': 'chr1', 'start': 1, 'end': 3, 'strand': '+', 'feature_type': 'promoter', 'coordinate_system': '1-based-closed', 'upstream': 2, 'downstream': 1}], 'promoter_count'...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 15) merge_gene_models

- **Purpose:** Document a reusable merge gene models capability for Annotation workflows. Candidate CAP-0049 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge gene models operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/annotation/local.py` and is exported through the formal DNAKit namespace. Coordinate conventions, contig lengths, and merge rules are explicit options; no annotation database or reference genome is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.merge_gene_models(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `merge_gene_models` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit coordinate provenance. It performs deterministic interval arithmetic only; no annotation database or reference sequence is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import merge_gene_models

result = merge_gene_models(
    input={"segments": [{"gene_id": "g1", "chrom": "chr1", "start": 1, "end": 3}, {"gene_id": "g1", "chrom": "chr1", "start": 4, "end": 5}]} ,
    options={"coordinate_system": "1-based-closed", "merge_adjacent": True}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0049",
  "operation": "merge_gene_models",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'merge_gene_models', 'result': {'models': [{'gene_id': 'g1', 'gene_name': 'g1', 'chrom': 'chr1', 'strand': '+', 'start': 1, 'end': 5, 'segments': [{'start': 1, 'end': 5}], 'segment_count': 1, 'coordinate_system': '1-based-closed'}], 'model_count': 1, 'input_se...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 16) parse

- **Purpose:** Document a reusable parse capability for Annotation workflows. Candidate CAP-0050 comes from 28 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the parse operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (cds, centromere, enhancer, exon, gene, gene_family, insulator, intron, lncrna, mirna, orf, ortholog, and 16 more). The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.parse(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `parse` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import parse

result = parse(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0050",
  "operation": "parse",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'parse', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed4389481...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 17) summarize

- **Purpose:** Document a reusable summarize capability for Annotation workflows. Candidate CAP-0051 comes from 28 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the summarize operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (cds, centromere, enhancer, exon, gene, gene_family, insulator, intron, lncrna, mirna, orf, ortholog, and 16 more). The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.annotation.summarize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `summarize` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import summarize

result = summarize(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0051",
  "operation": "summarize",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'summarize', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed438...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 18) validate

- **Purpose:** Document a reusable validate capability for Annotation workflows. Candidate CAP-0052 comes from 28 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (cds, centromere, enhancer, exon, gene, gene_family, insulator, intron, lncrna, mirna, orf, ortholog, and 16 more). The implementation is available at `src/dnakit/annotation/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `generic_operation_name`. Existing-name evidence: `dnakit.standardize.validate|validate|validate_set`.

- **API:** `dnakit.annotation.validate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `validate` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.annotation import validate

result = validate(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0052",
  "operation": "validate",
  "result_type": "dict",
  "result_preview": "{'module': 'annotation', 'operation': 'validate', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed4389...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```
