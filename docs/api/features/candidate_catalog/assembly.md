# Assembly

This category contains **37** documented candidate capabilities; 37 have a locally callable formal entry (0 partial due to binary-format dependencies; 0 bounded external backends; 32 dependency-adapter only). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) add_genome_to_pangenome

- **Purpose:** Document a reusable add genome to pangenome capability for Assembly workflows. Candidate CAP-0053 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the add genome to pangenome operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.add_genome_to_pangenome(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `add_genome_to_pangenome` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import add_genome_to_pangenome

result = add_genome_to_pangenome(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0053",
  "operation": "add_genome_to_pangenome",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'add_genome_to_pangenome', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 2) build_chromosome_level_assembly

- **Purpose:** Document a reusable build chromosome level assembly capability for Assembly workflows. Candidate CAP-0054 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the build chromosome level assembly operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.build_chromosome_level_assembly(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `build_chromosome_level_assembly` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import build_chromosome_level_assembly

result = build_chromosome_level_assembly(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0054",
  "operation": "build_chromosome_level_assembly",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'build_chromosome_level_assembly', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 3) build_command

- **Purpose:** Document a reusable build command capability for Assembly workflows. Candidate CAP-0055 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the build command operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.build_command(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `build_command` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import build_command

result = build_command(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0055",
  "operation": "build_command",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'build_command', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed4...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 4) build_graph_genome

- **Purpose:** Document a reusable build graph genome capability for Assembly workflows. Candidate CAP-0056 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the build graph genome operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.build_graph_genome(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `build_graph_genome` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import build_graph_genome

result = build_graph_genome(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0056",
  "operation": "build_graph_genome",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'build_graph_genome', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 5) build_haplotype_assembly

- **Purpose:** Document a reusable build haplotype assembly capability for Assembly workflows. Candidate CAP-0057 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the build haplotype assembly operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.build_haplotype_assembly(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `build_haplotype_assembly` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import build_haplotype_assembly

result = build_haplotype_assembly(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0057",
  "operation": "build_haplotype_assembly",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'build_haplotype_assembly', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 6) build_pangenome

- **Purpose:** Document a reusable build pangenome capability for Assembly workflows. Candidate CAP-0058 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the build pangenome operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.build_pangenome(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `build_pangenome` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import build_pangenome

result = build_pangenome(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0058",
  "operation": "build_pangenome",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'build_pangenome', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfe...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 7) calculate_assembly_completeness

- **Purpose:** Document a reusable calculate assembly completeness capability for Assembly workflows. Candidate CAP-0059 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate assembly completeness operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.calculate_assembly_completeness(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_assembly_completeness` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_assembly_completeness

result = calculate_assembly_completeness(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0059",
  "operation": "calculate_assembly_completeness",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'calculate_assembly_completeness', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 8) calculate_busco_completeness

- **Purpose:** Document a reusable calculate busco completeness capability for Assembly workflows. Candidate CAP-0060 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate busco completeness operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.calculate_busco_completeness(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_busco_completeness` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_busco_completeness

result = calculate_busco_completeness(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0060",
  "operation": "calculate_busco_completeness",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'calculate_busco_completeness', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 9) calculate_l50

- **Purpose:** Document a reusable calculate l50 capability for Assembly workflows. Candidate CAP-0061 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate l50 operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/assembly/local.py` and is exported through the formal DNAKit namespace. It computes deterministic contiguity statistics from explicit contig lengths; no assembler, graph, reference genome, or BUSCO database is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.calculate_l50(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `calculate_l50` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit local provenance. The metric is computed from positive contig lengths only; no assembler, graph, reference genome, or BUSCO database is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_l50

result = calculate_l50(
    input={"lengths": [100, 60, 40]},
    options={}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0061",
  "operation": "calculate_l50",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_l50', 'result': {'l50': 1, 'n_at_percentile': 100, 'percentile': 50.0, 'contig_count': 3, 'total_length': 200}, 'provenance': {'implementation': 'dnakit.native_assembly_metrics', 'backend': 'local_dependency_free', 'scientific_boundary': 'contiguity...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 10) calculate_n50

- **Purpose:** Document a reusable calculate n50 capability for Assembly workflows. Candidate CAP-0062 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate n50 operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/assembly/local.py` and is exported through the formal DNAKit namespace. It computes deterministic contiguity statistics from explicit contig lengths; no assembler, graph, reference genome, or BUSCO database is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.calculate_n50(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `calculate_n50` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit local provenance. The metric is computed from positive contig lengths only; no assembler, graph, reference genome, or BUSCO database is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_n50

result = calculate_n50(
    input={"lengths": [100, 60, 40]},
    options={}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0062",
  "operation": "calculate_n50",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_n50', 'result': {'n50': 100, 'l50': 1, 'contig_count': 3, 'total_length': 200, 'percentile': 50.0, 'options_used': []}, 'provenance': {'implementation': 'dnakit.native_assembly_metrics', 'backend': 'local_dependency_free', 'scientific_boundary': 'co...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 11) calculate_ng50

- **Purpose:** Document a reusable calculate ng50 capability for Assembly workflows. Candidate CAP-0063 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate ng50 operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/assembly/local.py` and is exported through the formal DNAKit namespace. It computes deterministic contiguity statistics from explicit contig lengths; no assembler, graph, reference genome, or BUSCO database is queried. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.calculate_ng50(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `calculate_ng50` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit local provenance. The metric is computed from positive contig lengths only; no assembler, graph, reference genome, or BUSCO database is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_ng50

result = calculate_ng50(
    input={"lengths": [100, 60, 40]},
    options={"genome_size": 250}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0063",
  "operation": "calculate_ng50",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_ng50', 'result': {'ng': 60, 'ng50': 60, 'contigs_to_threshold': 2, 'percentile': 50.0, 'expected_genome_size': 250, 'assembly_length': 200}, 'provenance': {'implementation': 'dnakit.native_assembly_metrics', 'backend': 'local_dependency_free', 'scie...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 12) calculate_nx

- **Purpose:** Document a reusable calculate nx capability for Assembly workflows. Candidate CAP-0064 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate nx operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/assembly/local.py` and is exported through the formal DNAKit namespace. It computes deterministic contiguity statistics from explicit contig lengths; no assembler, graph, reference genome, or BUSCO database is queried. Review basis: `no_exact_public_match`. Existing-name evidence: `calculate_dn_ds`.

- **API:** `dnakit.assembly.calculate_nx(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `calculate_nx` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit local provenance. The metric is computed from positive contig lengths only; no assembler, graph, reference genome, or BUSCO database is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_nx

result = calculate_nx(
    input={"lengths": [100, 60, 40]},
    options={}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0064",
  "operation": "calculate_nx",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_nx', 'result': {'percentile': 50.0, 'nx': 100, 'contigs_to_threshold': 1, 'contig_count': 3, 'total_length': 200}, 'provenance': {'implementation': 'dnakit.native_assembly_metrics', 'backend': 'local_dependency_free', 'scientific_boundary': 'contigu...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 13) calculate_pangenome_openness

- **Purpose:** Document a reusable calculate pangenome openness capability for Assembly workflows. Candidate CAP-0065 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate pangenome openness operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.calculate_pangenome_openness(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_pangenome_openness` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_pangenome_openness

result = calculate_pangenome_openness(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0065",
  "operation": "calculate_pangenome_openness",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'calculate_pangenome_openness', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 14) calculate_stats

- **Purpose:** Document a reusable calculate stats capability for Assembly workflows. Candidate CAP-0066 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate stats operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The dependency-free implementation is available at `src/dnakit/assembly/local.py` and is exported through the formal DNAKit namespace. It computes deterministic contiguity statistics from explicit contig lengths; no assembler, graph, reference genome, or BUSCO database is queried. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.calculate_stats(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `calculate_stats` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit local provenance. The metric is computed from positive contig lengths only; no assembler, graph, reference genome, or BUSCO database is queried.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import calculate_stats

result = calculate_stats(
    input={"lengths": [100, 60, 40]},
    options={}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0066",
  "operation": "calculate_stats",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_stats', 'result': {'contig_count': 3, 'total_length': 200, 'min_length': 40, 'max_length': 100, 'mean_length': 66.66666666666667, 'n50': 100, 'l50': 1, 'n90': 40, 'l90': 3, 'sorted_lengths': [100, 60, 40], 'options_used': []}, 'provenance': {'implem...",
  "provenance": {
    "implementation": "dnakit.native_assembly_metrics",
    "backend": "local_dependency_free",
    "scientific_boundary": "contiguity statistics over explicit contig lengths; no assembler or reference genome"
  }
}
```

## 15) check_input

- **Purpose:** Document a reusable check input capability for Assembly workflows. Candidate CAP-0067 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the check input operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.check_input(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `check_input` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import check_input

result = check_input(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0067",
  "operation": "check_input",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'check_input', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed438...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 16) compare_assemblies

- **Purpose:** Document a reusable compare assemblies capability for Assembly workflows. Candidate CAP-0068 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compare assemblies operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.compare_assemblies(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compare_assemblies` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import compare_assemblies

result = compare_assemblies(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0068",
  "operation": "compare_assemblies",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'compare_assemblies', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 17) compare_assembly_graphs

- **Purpose:** Document a reusable compare assembly graphs capability for Assembly workflows. Candidate CAP-0069 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compare assembly graphs operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.compare_assembly_graphs(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compare_assembly_graphs` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import compare_assembly_graphs

result = compare_assembly_graphs(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0069",
  "operation": "compare_assembly_graphs",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'compare_assembly_graphs', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 18) detect_assembly_contamination

- **Purpose:** Document a reusable detect assembly contamination capability for Assembly workflows. Candidate CAP-0070 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect assembly contamination operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.detect_assembly_contamination(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `detect_assembly_contamination` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import detect_assembly_contamination

result = detect_assembly_contamination(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0070",
  "operation": "detect_assembly_contamination",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'detect_assembly_contamination', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 19) detect_chimeric_contigs

- **Purpose:** Document a reusable detect chimeric contigs capability for Assembly workflows. Candidate CAP-0071 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect chimeric contigs operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.detect_chimeric_contigs(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `detect_chimeric_contigs` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import detect_chimeric_contigs

result = detect_chimeric_contigs(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0071",
  "operation": "detect_chimeric_contigs",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'detect_chimeric_contigs', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 20) detect_misassemblies

- **Purpose:** Document a reusable detect misassemblies capability for Assembly workflows. Candidate CAP-0072 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect misassemblies operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.detect_misassemblies(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `detect_misassemblies` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import detect_misassemblies

result = detect_misassemblies(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0072",
  "operation": "detect_misassemblies",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'detect_misassemblies', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 21) estimate_heterozygosity

- **Purpose:** Document a reusable estimate heterozygosity capability for Assembly workflows. Candidate CAP-0073 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the estimate heterozygosity operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.estimate_heterozygosity(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `estimate_heterozygosity` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import estimate_heterozygosity

result = estimate_heterozygosity(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0073",
  "operation": "estimate_heterozygosity",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'estimate_heterozygosity', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 22) estimate_repeats

- **Purpose:** Document a reusable estimate repeats capability for Assembly workflows. Candidate CAP-0074 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the estimate repeats operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.estimate_repeats(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `estimate_repeats` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import estimate_repeats

result = estimate_repeats(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0074",
  "operation": "estimate_repeats",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'estimate_repeats', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bf...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 23) estimate_size

- **Purpose:** Document a reusable estimate size capability for Assembly workflows. Candidate CAP-0075 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the estimate size operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.estimate_size(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `estimate_size` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import estimate_size

result = estimate_size(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0075",
  "operation": "estimate_size",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'estimate_size', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed4...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 24) extract_accessory_pangenome

- **Purpose:** Document a reusable extract accessory pangenome capability for Assembly workflows. Candidate CAP-0076 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract accessory pangenome operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.extract_accessory_pangenome(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `extract_accessory_pangenome` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import extract_accessory_pangenome

result = extract_accessory_pangenome(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0076",
  "operation": "extract_accessory_pangenome",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'extract_accessory_pangenome', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 25) extract_core_pangenome

- **Purpose:** Document a reusable extract core pangenome capability for Assembly workflows. Candidate CAP-0077 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract core pangenome operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.extract_core_pangenome(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `extract_core_pangenome` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import extract_core_pangenome

result = extract_core_pangenome(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0077",
  "operation": "extract_core_pangenome",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'extract_core_pangenome', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f8...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 26) identify_centromeres

- **Purpose:** Document a reusable identify centromeres capability for Assembly workflows. Candidate CAP-0078 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the identify centromeres operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.identify_centromeres(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `identify_centromeres` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import identify_centromeres

result = identify_centromeres(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0078",
  "operation": "identify_centromeres",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'identify_centromeres', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 27) identify_repeats_in_assembly

- **Purpose:** Document a reusable identify repeats in assembly capability for Assembly workflows. Candidate CAP-0079 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the identify repeats in assembly operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.identify_repeats_in_assembly(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `identify_repeats_in_assembly` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import identify_repeats_in_assembly

result = identify_repeats_in_assembly(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0079",
  "operation": "identify_repeats_in_assembly",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'identify_repeats_in_assembly', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 28) identify_telomeres

- **Purpose:** Document a reusable identify telomeres capability for Assembly workflows. Candidate CAP-0080 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the identify telomeres operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.identify_telomeres(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `identify_telomeres` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import identify_telomeres

result = identify_telomeres(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0080",
  "operation": "identify_telomeres",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'identify_telomeres', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 29) parse_log

- **Purpose:** Document a reusable parse log capability for Assembly workflows. Candidate CAP-0081 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the parse log operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.parse_log(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `parse_log` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import parse_log

result = parse_log(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0081",
  "operation": "parse_log",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'parse_log', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 30) parse_output

- **Purpose:** Document a reusable parse output capability for Assembly workflows. Candidate CAP-0082 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the parse output operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.parse_output(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `parse_output` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import parse_output

result = parse_output(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0082",
  "operation": "parse_output",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'parse_output', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 31) phase_assembly

- **Purpose:** Document a reusable phase assembly capability for Assembly workflows. Candidate CAP-0083 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the phase assembly operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.phase_assembly(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `phase_assembly` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import phase_assembly

result = phase_assembly(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0083",
  "operation": "phase_assembly",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'phase_assembly', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 32) plot_graph

- **Purpose:** Document a reusable plot graph capability for Assembly workflows. Candidate CAP-0084 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the plot graph operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.plot_graph(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `plot_graph` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import plot_graph

result = plot_graph(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0084",
  "operation": "plot_graph",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'plot_graph', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed4389...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 33) polish

- **Purpose:** Document a reusable polish capability for Assembly workflows. Candidate CAP-0085 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the polish operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.polish(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `polish` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import polish

result = polish(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0085",
  "operation": "polish",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'polish', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894816...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 34) purge_haplotigs

- **Purpose:** Document a reusable purge haplotigs capability for Assembly workflows. Candidate CAP-0086 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the purge haplotigs operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.purge_haplotigs(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `purge_haplotigs` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import purge_haplotigs

result = purge_haplotigs(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0086",
  "operation": "purge_haplotigs",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'purge_haplotigs', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfe...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 35) remove_contamination

- **Purpose:** Document a reusable remove contamination capability for Assembly workflows. Candidate CAP-0087 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the remove contamination operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.remove_contamination(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `remove_contamination` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import remove_contamination

result = remove_contamination(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0087",
  "operation": "remove_contamination",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'remove_contamination', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 36) run

- **Purpose:** Document a reusable run capability for Assembly workflows. Candidate CAP-0088 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the run operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `generic_operation_name`. Existing-name evidence: `dnakit.workflows.run|run`.

- **API:** `dnakit.assembly.run(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `run` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import run

result = run(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0088",
  "operation": "run",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'run', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894816a'}...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 37) scaffold

- **Purpose:** Document a reusable scaffold capability for Assembly workflows. Candidate CAP-0089 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the scaffold operation, and return structured results with provenance and diagnostics. The optional `assembler` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/assembly/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.assembly.scaffold(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `scaffold` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.assembly import scaffold

result = scaffold(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0089",
  "operation": "scaffold",
  "result_type": "dict",
  "result_preview": "{'module': 'assembly', 'operation': 'scaffold', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed438948...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```
