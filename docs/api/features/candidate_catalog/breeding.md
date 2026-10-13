# Breeding

This category contains **7** documented candidate capabilities; 7 have a locally callable formal entry (0 partial due to binary-format dependencies; 0 bounded external backends; 7 dependency-adapter only). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) calculate

- **Purpose:** Document a reusable calculate capability for Breeding workflows. Candidate CAP-0090 comes from 15 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (conservation, disease_resistance, gblup, gebv, genomic_relationship, genotype_environment, haplotype_breeding, heterosis, inbreeding, marker_assisted_selection, mating, multi_trait, and 3 more). The implementation is available at `src/dnakit/breeding/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.breeding.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.breeding import calculate

result = calculate(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0090",
  "operation": "calculate",
  "result_type": "dict",
  "result_preview": "{'module': 'breeding', 'operation': 'calculate', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 2) fit

- **Purpose:** Document a reusable fit capability for Breeding workflows. Candidate CAP-0092 comes from 17 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the fit operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (conservation, disease_resistance, gblup, gebv, genomic_relationship, genotype_environment, gwas, haplotype_breeding, heterosis, inbreeding, marker_assisted_selection, mating, and 5 more). The implementation is available at `src/dnakit/breeding/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.breeding.fit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `fit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.breeding import fit

result = fit(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0092",
  "operation": "fit",
  "result_type": "dict",
  "result_preview": "{'module': 'breeding', 'operation': 'fit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894816a'}...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 3) optimize

- **Purpose:** Document a reusable optimize capability for Breeding workflows. Candidate CAP-0093 comes from 17 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the optimize operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (conservation, disease_resistance, gblup, gebv, genomic_relationship, genotype_environment, gwas, haplotype_breeding, heterosis, inbreeding, marker_assisted_selection, mating, and 5 more). The implementation is available at `src/dnakit/breeding/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.breeding.optimize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `optimize` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.breeding import optimize

result = optimize(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0093",
  "operation": "optimize",
  "result_type": "dict",
  "result_preview": "{'module': 'breeding', 'operation': 'optimize', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed438948...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 4) plot

- **Purpose:** Document a reusable plot capability for Breeding workflows. Candidate CAP-0094 comes from 15 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the plot operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (conservation, disease_resistance, gblup, gebv, genomic_relationship, genotype_environment, haplotype_breeding, heterosis, inbreeding, marker_assisted_selection, mating, multi_trait, and 3 more). The implementation is available at `src/dnakit/breeding/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.breeding.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `plot` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.breeding import plot

result = plot(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0094",
  "operation": "plot",
  "result_type": "dict",
  "result_preview": "{'module': 'breeding', 'operation': 'plot', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894816a'...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 5) predict

- **Purpose:** Document a reusable predict capability for Breeding workflows. Candidate CAP-0095 comes from 17 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the predict operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (conservation, disease_resistance, gblup, gebv, genomic_relationship, genotype_environment, gwas, haplotype_breeding, heterosis, inbreeding, marker_assisted_selection, mating, and 5 more). The implementation is available at `src/dnakit/breeding/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.breeding.predict(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `predict` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.breeding import predict

result = predict(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0095",
  "operation": "predict",
  "result_type": "dict",
  "result_preview": "{'module': 'breeding', 'operation': 'predict', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed4389481...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 6) select

- **Purpose:** Document a reusable select capability for Breeding workflows. Candidate CAP-0096 comes from 17 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the select operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (conservation, disease_resistance, gblup, gebv, genomic_relationship, genotype_environment, gwas, haplotype_breeding, heterosis, inbreeding, marker_assisted_selection, mating, and 5 more). The implementation is available at `src/dnakit/breeding/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.breeding.select(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `select` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.breeding import select

result = select(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0096",
  "operation": "select",
  "result_type": "dict",
  "result_preview": "{'module': 'breeding', 'operation': 'select', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed43894816...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 7) simulate

- **Purpose:** Document a reusable simulate capability for Breeding workflows. Candidate CAP-0097 comes from 17 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the simulate operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (conservation, disease_resistance, gblup, gebv, genomic_relationship, genotype_environment, gwas, haplotype_breeding, heterosis, inbreeding, marker_assisted_selection, mating, and 5 more). The implementation is available at `src/dnakit/breeding/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. Existing-name evidence: `simulate_pcr`.

- **API:** `dnakit.breeding.simulate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `simulate` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.breeding import simulate

result = simulate(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0097",
  "operation": "simulate",
  "result_type": "dict",
  "result_preview": "{'module': 'breeding', 'operation': 'simulate', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed438948...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```
