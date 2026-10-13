# Replication and DNA repair

This category contains **10** documented candidate capabilities; 10 have a locally callable formal entry (0 partial due to binary-format dependencies; 0 bounded external backends; 10 dependency-adapter only). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) annotate

- **Purpose:** Document a reusable annotate capability for Replication and DNA repair workflows. Candidate CAP-0573 comes from 15 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the annotate operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (crossover, damage_hotspot, double_strand_break, gene_conversion, microhomology, mutational_signature, okazaki_fragment, recombination_hotspot, repair_breakpoint, repair_genes, replication_fork, replication_origin, and 3 more). The implementation is available at `src/dnakit/replication_repair/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.replication_repair.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `annotate` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.replication_repair import annotate

result = annotate(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0573",
  "operation": "annotate",
  "result_type": "dict",
  "result_preview": "{'module': 'replication_repair', 'operation': 'annotate', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 2) build_recombination_map

- **Purpose:** Document a reusable build recombination map capability for Replication and DNA repair workflows. Candidate CAP-0574 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the build recombination map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/replication_repair/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.replication_repair.build_recombination_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `build_recombination_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.replication_repair import build_recombination_map

result = build_recombination_map(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0574",
  "operation": "build_recombination_map",
  "result_type": "dict",
  "result_preview": "{'module': 'replication_repair', 'operation': 'build_recombination_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49d...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 3) calculate

- **Purpose:** Document a reusable calculate capability for Replication and DNA repair workflows. Candidate CAP-0575 comes from 18 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (centromere, crossover, damage_hotspot, double_strand_break, gene_conversion, microhomology, mutational_signature, mutational_signature_exposure, okazaki_fragment, recombination_hotspot, recombination_rate, repair_breakpoint, and 6 more). The implementation is available at `src/dnakit/replication_repair/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.replication_repair.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.replication_repair import calculate

result = calculate(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0575",
  "operation": "calculate",
  "result_type": "dict",
  "result_preview": "{'module': 'replication_repair', 'operation': 'calculate', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 4) classify_repair_pathway

- **Purpose:** Document a reusable classify repair pathway capability for Replication and DNA repair workflows. Candidate CAP-0576 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the classify repair pathway operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/replication_repair/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.replication_repair.classify_repair_pathway(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `classify_repair_pathway` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.replication_repair import classify_repair_pathway

result = classify_repair_pathway(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0576",
  "operation": "classify_repair_pathway",
  "result_type": "dict",
  "result_preview": "{'module': 'replication_repair', 'operation': 'classify_repair_pathway', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49d...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 5) compare

- **Purpose:** Document a reusable compare capability for Replication and DNA repair workflows. Candidate CAP-0577 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compare operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (crossover, damage_hotspot, double_strand_break, gene_conversion, microhomology, mutational_signature, okazaki_fragment, recombination_hotspot, repair_breakpoint, replication_fork, replication_origin, replication_stress, and 2 more). The implementation is available at `src/dnakit/replication_repair/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `generic_operation_name`. Existing-name evidence: `compare|dnakit.similarity.compare`.

- **API:** `dnakit.replication_repair.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compare` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.replication_repair import compare

result = compare(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0577",
  "operation": "compare",
  "result_type": "dict",
  "result_preview": "{'module': 'replication_repair', 'operation': 'compare', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 6) detect

- **Purpose:** Document a reusable detect capability for Replication and DNA repair workflows. Candidate CAP-0578 comes from 16 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (centromere, crossover, damage_hotspot, double_strand_break, gene_conversion, microhomology, mutational_signature, okazaki_fragment, recombination_hotspot, repair_breakpoint, replication_fork, replication_origin, and 4 more). The implementation is available at `src/dnakit/replication_repair/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.replication_repair.detect(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `detect` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.replication_repair import detect

result = detect(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0578",
  "operation": "detect",
  "result_type": "dict",
  "result_preview": "{'module': 'replication_repair', 'operation': 'detect', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bf...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 7) map

- **Purpose:** Document a reusable map capability for Replication and DNA repair workflows. Candidate CAP-0579 comes from 16 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the map operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (centromere, crossover, damage_hotspot, double_strand_break, gene_conversion, microhomology, mutational_signature, okazaki_fragment, recombination_hotspot, repair_breakpoint, replication_fork, replication_origin, and 4 more). The implementation is available at `src/dnakit/replication_repair/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.replication_repair.map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.replication_repair import map

result = map(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0579",
  "operation": "map",
  "result_type": "dict",
  "result_preview": "{'module': 'replication_repair', 'operation': 'map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed4...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 8) plot

- **Purpose:** Document a reusable plot capability for Replication and DNA repair workflows. Candidate CAP-0580 comes from 16 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the plot operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (centromere, crossover, damage_hotspot, double_strand_break, gene_conversion, microhomology, mutational_signature, okazaki_fragment, recombination_hotspot, repair_breakpoint, replication_fork, replication_origin, and 4 more). The implementation is available at `src/dnakit/replication_repair/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.replication_repair.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `plot` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.replication_repair import plot

result = plot(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0580",
  "operation": "plot",
  "result_type": "dict",
  "result_preview": "{'module': 'replication_repair', 'operation': 'plot', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6bfed...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 9) predict

- **Purpose:** Document a reusable predict capability for Replication and DNA repair workflows. Candidate CAP-0581 comes from 19 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the predict operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (centromere, crossover, damage_hotspot, double_strand_break, gene_conversion, hr_repair_outcome, microhomology, mmej_repair_outcome, mutational_signature, nhej_repair_outcome, okazaki_fragment, recombination_hotspot, and 7 more). The implementation is available at `src/dnakit/replication_repair/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.replication_repair.predict(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `predict` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.replication_repair import predict

result = predict(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0581",
  "operation": "predict",
  "result_type": "dict",
  "result_preview": "{'module': 'replication_repair', 'operation': 'predict', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 10) simulate

- **Purpose:** Document a reusable simulate capability for Replication and DNA repair workflows. Candidate CAP-0582 comes from 16 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the simulate operation, and return structured results with provenance and diagnostics. The optional `entity` parameter selects one supported variant (centromere, crossover, damage_hotspot, double_strand_break, gene_conversion, microhomology, mutational_signature, okazaki_fragment, recombination_hotspot, repair_breakpoint, replication_fork, replication_origin, and 4 more). The implementation is available at `src/dnakit/replication_repair/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. Existing-name evidence: `simulate_pcr`.

- **API:** `dnakit.replication_repair.simulate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `simulate` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.replication_repair import simulate

result = simulate(
    input='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0582",
  "operation": "simulate",
  "result_type": "dict",
  "result_preview": "{'module': 'replication_repair', 'operation': 'simulate', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 4, 'text_length': 16, 'sequence_bases': 8, 'lengths': [4, 4], 'sha256': 'b84d5f88bf0b633a597cb3b696489b49a49df7cdd6021f872a6...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```
