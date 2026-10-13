# FASTA and FASTQ

This category contains **66** documented candidate capabilities; 66 have a locally callable formal entry (0 partial due to binary-format dependencies; 0 bounded external backends; 0 dependency-adapter only). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) build_fasta_index

- **Purpose:** Document a reusable build fasta index capability for FASTA and FASTQ workflows. Candidate CAP-0176 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the build fasta index operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `generic_operation_name`. Existing-name evidence: `build_fasta_index|build_fastq_index|build_sketch_index|dnakit.io.build_fasta_index|load_fasta_index`.

- **API:** `dnakit.fasta_fastq.build_fasta_index(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `build_fasta_index` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import build_fasta_index

result = build_fasta_index(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0176",
  "operation": "build_fasta_index",
  "result_type": "dict",
  "result_preview": "{'a': {'record_index': 0, 'length': 4}, 'b': {'record_index': 1, 'length': 4}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 2) build_fastq_index

- **Purpose:** Document a reusable build fastq index capability for FASTA and FASTQ workflows. Candidate CAP-0177 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the build fastq index operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `generic_operation_name`. Existing-name evidence: `build_fasta_index|build_fastq_index|build_sketch_index|dnakit.io.build_fastq_index|load_fastq_index`.

- **API:** `dnakit.fasta_fastq.build_fastq_index(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `build_fastq_index` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import build_fastq_index

result = build_fastq_index(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0177",
  "operation": "build_fastq_index",
  "result_type": "dict",
  "result_preview": "{'a': {'record_index': 0, 'length': 4, 'quality_length': 4}, 'b': {'record_index': 1, 'length': 4, 'quality_length': 4}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 3) calculate_gc_profile

- **Purpose:** Document a reusable calculate gc profile capability for FASTA and FASTQ workflows. Candidate CAP-0178 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate gc profile operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.calculate_gc_profile(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `calculate_gc_profile` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import calculate_gc_profile

result = calculate_gc_profile(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0178",
  "operation": "calculate_gc_profile",
  "result_type": "list",
  "result_preview": "[{'id': 'a', 'length': 4, 'gc_fraction': 0.5, 'gc_percent': 50.0}, {'id': 'b', 'length': 4, 'gc_fraction': 0.5, 'gc_percent': 50.0}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 4) calculate_length_distribution

- **Purpose:** Document a reusable calculate length distribution capability for FASTA and FASTQ workflows. Candidate CAP-0179 comes from 2 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate length distribution operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta, fastq). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.calculate_length_distribution(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `calculate_length_distribution` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import calculate_length_distribution

result = calculate_length_distribution(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0179",
  "operation": "calculate_length_distribution",
  "result_type": "dict",
  "result_preview": "{'count': 2, 'min': 4, 'max': 4, 'mean': 4.0, 'counts': {4: 2}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 5) calculate_n_content_profile

- **Purpose:** Document a reusable calculate n content profile capability for FASTA and FASTQ workflows. Candidate CAP-0180 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate n content profile operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.calculate_n_content_profile(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `calculate_n_content_profile` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import calculate_n_content_profile

result = calculate_n_content_profile(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0180",
  "operation": "calculate_n_content_profile",
  "result_type": "list",
  "result_preview": "[{'id': 'a', 'n_count': 0, 'n_fraction': 0.0}, {'id': 'b', 'n_count': 0, 'n_fraction': 0.0}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 6) calculate_quality_distribution

- **Purpose:** Document a reusable calculate quality distribution capability for FASTA and FASTQ workflows. Candidate CAP-0181 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate quality distribution operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.calculate_quality_distribution(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `calculate_quality_distribution` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import calculate_quality_distribution

result = calculate_quality_distribution(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0181",
  "operation": "calculate_quality_distribution",
  "result_type": "dict",
  "result_preview": "{'count': 8, 'min': 40, 'max': 40, 'mean': 40.0, 'counts': {40: 8}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 7) compare_fasta_files

- **Purpose:** Document a reusable compare fasta files capability for FASTA and FASTQ workflows. Candidate CAP-0182 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compare fasta files operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.compare_fasta_files(left[required], right[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `left`, `right`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `compare_fasta_files` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import compare_fasta_files

result = compare_fasta_files(
    left='>a\nACGT\n>b\nGGAA\n',
    right='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0182",
  "operation": "compare_fasta_files",
  "result_type": "dict",
  "result_preview": "{'left_only': [], 'right_only': [], 'common': ['a', 'b'], 'different_sequences': [], 'identical': True}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 8) convert_from_fasta

- **Purpose:** Document a reusable convert from fasta capability for FASTA and FASTQ workflows. Candidate CAP-0183 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert from fasta operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.convert_from_fasta(source[required], quality[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: `quality`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `convert_from_fasta` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import convert_from_fasta

result = convert_from_fasta(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0183",
  "operation": "convert_from_fasta",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 9) convert_from_tabular

- **Purpose:** Document a reusable convert from tabular capability for FASTA and FASTQ workflows. Candidate CAP-0184 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert from tabular operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.convert_from_tabular(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `convert_from_tabular` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import convert_from_tabular

result = convert_from_tabular(
    source='a\tACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0184",
  "operation": "convert_from_tabular",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 10) convert_to_fasta

- **Purpose:** Document a reusable convert to fasta capability for FASTA and FASTQ workflows. Candidate CAP-0185 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert to fasta operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.convert_to_fasta(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `convert_to_fasta` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import convert_to_fasta

result = convert_to_fasta(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0185",
  "operation": "convert_to_fasta",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 11) convert_to_tabular

- **Purpose:** Document a reusable convert to tabular capability for FASTA and FASTQ workflows. Candidate CAP-0186 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert to tabular operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.convert_to_tabular(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `convert_to_tabular` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import convert_to_tabular

result = convert_to_tabular(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0186",
  "operation": "convert_to_tabular",
  "result_type": "list",
  "result_preview": "[{'id': 'a', 'sequence': 'ACGT'}, {'id': 'b', 'sequence': 'GGAA'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 12) count_reads

- **Purpose:** Document a reusable count reads capability for FASTA and FASTQ workflows. Candidate CAP-0187 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the count reads operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.count_reads(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `count_reads` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import count_reads

result = count_reads(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0187",
  "operation": "count_reads",
  "result_type": "int",
  "result_preview": "2",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 13) count_records

- **Purpose:** Document a reusable count records capability for FASTA and FASTQ workflows. Candidate CAP-0188 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the count records operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.count_records(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `count_records` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import count_records

result = count_records(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0188",
  "operation": "count_records",
  "result_type": "int",
  "result_preview": "2",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 14) deinterleave_pairs

- **Purpose:** Document a reusable deinterleave pairs capability for FASTA and FASTQ workflows. Candidate CAP-0189 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the deinterleave pairs operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.deinterleave_pairs(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `deinterleave_pairs` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import deinterleave_pairs

result = deinterleave_pairs(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0189",
  "operation": "deinterleave_pairs",
  "result_type": "tuple",
  "result_preview": "([FastqRecord(identifier='a', sequence='ACGT', quality='IIII')], [FastqRecord(identifier='b', sequence='GGAA', quality='IIII')])",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 15) estimate_duplication

- **Purpose:** Document a reusable estimate duplication capability for FASTA and FASTQ workflows. Candidate CAP-0190 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the estimate duplication operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.estimate_duplication(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `estimate_duplication` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import estimate_duplication

result = estimate_duplication(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0190",
  "operation": "estimate_duplication",
  "result_type": "dict",
  "result_preview": "{'records': 2, 'unique_sequences': 2, 'duplicate_records': 0, 'duplicate_fraction': 0.0}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 16) estimate_error_rate

- **Purpose:** Document a reusable estimate error rate capability for FASTA and FASTQ workflows. Candidate CAP-0191 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the estimate error rate operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.estimate_error_rate(source[required], expected_quality[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: `expected_quality`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `estimate_error_rate` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import estimate_error_rate

result = estimate_error_rate(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0191",
  "operation": "estimate_error_rate",
  "result_type": "float",
  "result_preview": "0.0",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 17) extract_cds

- **Purpose:** Document a reusable extract cds capability for FASTA and FASTQ workflows. Candidate CAP-0192 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract cds operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.extract_cds(source[required], identifier[required], intervals[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `identifier`, `intervals`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `extract_cds` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_cds

result = extract_cds(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a",
    intervals=[(0, 3)]
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0192",
  "operation": "extract_cds",
  "result_type": "tuple",
  "result_preview": "('a:cds', 'ACG')",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 18) extract_exon

- **Purpose:** Document a reusable extract exon capability for FASTA and FASTQ workflows. Candidate CAP-0193 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract exon operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.extract_exon(source[required], identifier[required], start[required], end[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `identifier`, `start`, `end`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `extract_exon` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_exon

result = extract_exon(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a",
    start=0,
    end=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0193",
  "operation": "extract_exon",
  "result_type": "tuple",
  "result_preview": "('a:0-3', 'ACG')",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 19) extract_flanking

- **Purpose:** Document a reusable extract flanking capability for FASTA and FASTQ workflows. Candidate CAP-0194 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract flanking operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.extract_flanking(source[required], identifier[required], start[required], end[required], flank[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `identifier`, `start`, `end`; optional arguments: `flank`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `extract_flanking` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_flanking

result = extract_flanking(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a",
    start=0,
    end=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0194",
  "operation": "extract_flanking",
  "result_type": "tuple",
  "result_preview": "('a:flank', 'ACGT')",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 20) extract_introns

- **Purpose:** Document a reusable extract introns capability for FASTA and FASTQ workflows. Candidate CAP-0195 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract introns operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.extract_introns(source[required], identifier[required], exons[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `identifier`, `exons`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `extract_introns` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_introns

result = extract_introns(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a",
    exons=[(0, 2), (4, 6)]
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0195",
  "operation": "extract_introns",
  "result_type": "list",
  "result_preview": "[('a:intron:1', 'GT')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 21) extract_promoter

- **Purpose:** Document a reusable extract promoter capability for FASTA and FASTQ workflows. Candidate CAP-0196 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract promoter operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.extract_promoter(source[required], identifier[required], tss[required], upstream[optional], downstream[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `identifier`, `tss`; optional arguments: `upstream`, `downstream`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `extract_promoter` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_promoter

result = extract_promoter(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a",
    tss=2
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0196",
  "operation": "extract_promoter",
  "result_type": "tuple",
  "result_preview": "('a:promoter', 'ACGT')",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 22) extract_read_names

- **Purpose:** Document a reusable extract read names capability for FASTA and FASTQ workflows. Candidate CAP-0197 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract read names operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.extract_read_names(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `extract_read_names` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_read_names

result = extract_read_names(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0197",
  "operation": "extract_read_names",
  "result_type": "list",
  "result_preview": "['a', 'b']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 23) extract_region

- **Purpose:** Document a reusable extract region capability for FASTA and FASTQ workflows. Candidate CAP-0198 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract region operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.extract_region(source[required], identifier[required], start[required], end[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `identifier`, `start`, `end`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `extract_region` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_region

result = extract_region(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a",
    start=0,
    end=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0198",
  "operation": "extract_region",
  "result_type": "tuple",
  "result_preview": "('a:0-3', 'ACG')",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 24) extract_regions_from_bed

- **Purpose:** Document a reusable extract regions from bed capability for FASTA and FASTQ workflows. Candidate CAP-0199 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the extract regions from bed operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.extract_regions_from_bed(source[required], bed[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `bed`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `extract_regions_from_bed` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import extract_regions_from_bed

result = extract_regions_from_bed(
    source='>a\nACGT\n>b\nGGAA\n',
    bed='a\t0\t3\tregion1\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0199",
  "operation": "extract_regions_from_bed",
  "result_type": "list",
  "result_preview": "[('region1', 'ACG')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 25) filter

- **Purpose:** Document a reusable filter capability for FASTA and FASTQ workflows. Candidate CAP-0200 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.filter(source[required], min_length[optional], max_length[optional], ids[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: `min_length`, `max_length`, `ids`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `filter` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import filter

result = filter(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0200",
  "operation": "filter",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 26) filter_by_length

- **Purpose:** Document a reusable filter by length capability for FASTA and FASTQ workflows. Candidate CAP-0201 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter by length operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.filter_by_length(source[required], min_length[optional], max_length[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: `min_length`, `max_length`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `filter_by_length` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import filter_by_length

result = filter_by_length(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0201",
  "operation": "filter_by_length",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 27) filter_by_quality

- **Purpose:** Document a reusable filter by quality capability for FASTA and FASTQ workflows. Candidate CAP-0202 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter by quality operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.filter_by_quality(source[required], min_quality[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `min_quality`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `filter_by_quality` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import filter_by_quality

result = filter_by_quality(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    min_quality=20
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0202",
  "operation": "filter_by_quality",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 28) get_ids

- **Purpose:** Document a reusable get ids capability for FASTA and FASTQ workflows. Candidate CAP-0203 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the get ids operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.get_ids(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `get_ids` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import get_ids

result = get_ids(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0203",
  "operation": "get_ids",
  "result_type": "list",
  "result_preview": "['a', 'b']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 29) get_lengths

- **Purpose:** Document a reusable get lengths capability for FASTA and FASTQ workflows. Candidate CAP-0204 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the get lengths operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.get_lengths(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `get_lengths` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import get_lengths

result = get_lengths(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0204",
  "operation": "get_lengths",
  "result_type": "dict",
  "result_preview": "{'a': 4, 'b': 4}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 30) get_sequence_by_id

- **Purpose:** Document a reusable get sequence by id capability for FASTA and FASTQ workflows. Candidate CAP-0205 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the get sequence by id operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. Existing-name evidence: `sequence_by_id`.

- **API:** `dnakit.fasta_fastq.get_sequence_by_id(source[required], identifier[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `identifier`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `get_sequence_by_id` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import get_sequence_by_id

result = get_sequence_by_id(
    source='>a\nACGT\n>b\nGGAA\n',
    identifier="a"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0205",
  "operation": "get_sequence_by_id",
  "result_type": "str",
  "result_preview": "'ACGT'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 31) interleave_pairs

- **Purpose:** Document a reusable interleave pairs capability for FASTA and FASTQ workflows. Candidate CAP-0206 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the interleave pairs operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.interleave_pairs(first[required], second[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `first`, `second`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `interleave_pairs` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import interleave_pairs

result = interleave_pairs(
    first='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    second='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0206",
  "operation": "interleave_pairs",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 32) iterate

- **Purpose:** Document a reusable iterate capability for FASTA and FASTQ workflows. Candidate CAP-0207 comes from 2 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the iterate operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta, fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. Existing-name evidence: `literature`.

- **API:** `dnakit.fasta_fastq.iterate(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `iterate` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import iterate

result = iterate(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0207",
  "operation": "iterate",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 33) mask_regions

- **Purpose:** Document a reusable mask regions capability for FASTA and FASTQ workflows. Candidate CAP-0208 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the mask regions operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.mask_regions(source[required], regions[required], symbol[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `regions`; optional arguments: `symbol`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `mask_regions` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import mask_regions

result = mask_regions(
    source='>a\nACGT\n>b\nGGAA\n',
    regions={"a": [(0, 2)]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0208",
  "operation": "mask_regions",
  "result_type": "list",
  "result_preview": "[('a', 'NNGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 34) merge

- **Purpose:** Document a reusable merge capability for FASTA and FASTQ workflows. Candidate CAP-0209 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.merge(sources[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sources`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `merge` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import merge

result = merge(
    sources=['>a\nACGT\n>b\nGGAA\n']
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0209",
  "operation": "merge",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 35) merge_lanes

- **Purpose:** Document a reusable merge lanes capability for FASTA and FASTQ workflows. Candidate CAP-0210 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge lanes operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.merge_lanes(sources[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `sources`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `merge_lanes` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import merge_lanes

result = merge_lanes(
    sources=['@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n']
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0210",
  "operation": "merge_lanes",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 36) query_fasta_index

- **Purpose:** Document a reusable query fasta index capability for FASTA and FASTQ workflows. Candidate CAP-0211 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the query fasta index operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.query_fasta_index(index[required], identifier[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `index`, `identifier`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `query_fasta_index` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import query_fasta_index

result = query_fasta_index(
    index={"a": {"record_index": 0, "length": 4}},
    identifier="a"
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0211",
  "operation": "query_fasta_index",
  "result_type": "dict",
  "result_preview": "{'record_index': 0, 'length': 4}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 37) remove_duplicate_ids

- **Purpose:** Document a reusable remove duplicate ids capability for FASTA and FASTQ workflows. Candidate CAP-0212 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the remove duplicate ids operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.remove_duplicate_ids(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `remove_duplicate_ids` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import remove_duplicate_ids

result = remove_duplicate_ids(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0212",
  "operation": "remove_duplicate_ids",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 38) remove_duplicate_sequences

- **Purpose:** Document a reusable remove duplicate sequences capability for FASTA and FASTQ workflows. Candidate CAP-0213 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the remove duplicate sequences operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.remove_duplicate_sequences(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `remove_duplicate_sequences` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import remove_duplicate_sequences

result = remove_duplicate_sequences(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0213",
  "operation": "remove_duplicate_sequences",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 39) remove_poly_a

- **Purpose:** Document a reusable remove poly a capability for FASTA and FASTQ workflows. Candidate CAP-0214 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the remove poly a operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.remove_poly_a(source[required], minimum[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: `minimum`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `remove_poly_a` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import remove_poly_a

result = remove_poly_a(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0214",
  "operation": "remove_poly_a",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 40) remove_poly_g

- **Purpose:** Document a reusable remove poly g capability for FASTA and FASTQ workflows. Candidate CAP-0215 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the remove poly g operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.remove_poly_g(source[required], minimum[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: `minimum`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `remove_poly_g` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import remove_poly_g

result = remove_poly_g(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0215",
  "operation": "remove_poly_g",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 41) rename_headers

- **Purpose:** Document a reusable rename headers capability for FASTA and FASTQ workflows. Candidate CAP-0216 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the rename headers operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.rename_headers(source[required], prefix[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: `prefix`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `rename_headers` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import rename_headers

result = rename_headers(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0216",
  "operation": "rename_headers",
  "result_type": "list",
  "result_preview": "[('record_1', 'ACGT'), ('record_2', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 42) repair_fasta_headers

- **Purpose:** Document a reusable repair fasta headers capability for FASTA and FASTQ workflows. Candidate CAP-0217 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the repair fasta headers operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.repair_fasta_headers(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `repair_fasta_headers` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import repair_fasta_headers

result = repair_fasta_headers(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0217",
  "operation": "repair_fasta_headers",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 43) repair_pairs

- **Purpose:** Document a reusable repair pairs capability for FASTA and FASTQ workflows. Candidate CAP-0218 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the repair pairs operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.repair_pairs(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `repair_pairs` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import repair_pairs

result = repair_pairs(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0218",
  "operation": "repair_pairs",
  "result_type": "list",
  "result_preview": "[]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 44) reverse_complement_records

- **Purpose:** Document a reusable reverse complement records capability for FASTA and FASTQ workflows. Candidate CAP-0219 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the reverse complement records operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `merged_dimension`. Existing-name evidence: `reverse_complement|reverse_complement_record|reverse_complement_search`.

- **API:** `dnakit.fasta_fastq.reverse_complement_records(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `reverse_complement_records` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import reverse_complement_records

result = reverse_complement_records(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0219",
  "operation": "reverse_complement_records",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'TTCC')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 45) sample_records

- **Purpose:** Document a reusable sample records capability for FASTA and FASTQ workflows. Candidate CAP-0221 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sample records operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.sample_records(source[required], size[required], seed[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `size`; optional arguments: `seed`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `sample_records` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import sample_records

result = sample_records(
    source='>a\nACGT\n>b\nGGAA\n',
    size=1,
    seed=7
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0221",
  "operation": "sample_records",
  "result_type": "list",
  "result_preview": "[('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 46) sanitize_headers

- **Purpose:** Document a reusable sanitize headers capability for FASTA and FASTQ workflows. Candidate CAP-0222 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sanitize headers operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.sanitize_headers(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `sanitize_headers` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import sanitize_headers

result = sanitize_headers(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0222",
  "operation": "sanitize_headers",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 47) shuffle_reads

- **Purpose:** Document a reusable shuffle reads capability for FASTA and FASTQ workflows. Candidate CAP-0223 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the shuffle reads operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.shuffle_reads(source[required], seed[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: `seed`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `shuffle_reads` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import shuffle_reads

result = shuffle_reads(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    seed=7
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0223",
  "operation": "shuffle_reads",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 48) shuffle_records

- **Purpose:** Document a reusable shuffle records capability for FASTA and FASTQ workflows. Candidate CAP-0224 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the shuffle records operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.shuffle_records(source[required], seed[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: `seed`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `shuffle_records` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import shuffle_records

result = shuffle_records(
    source='>a\nACGT\n>b\nGGAA\n',
    seed=7
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0224",
  "operation": "shuffle_records",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 49) sort

- **Purpose:** Document a reusable sort capability for FASTA and FASTQ workflows. Candidate CAP-0225 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sort operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.sort(source[required], by[optional], reverse[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: `by`, `reverse`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `sort` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import sort

result = sort(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0225",
  "operation": "sort",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 50) split_by_barcode

- **Purpose:** Document a reusable split by barcode capability for FASTA and FASTQ workflows. Candidate CAP-0226 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split by barcode operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.split_by_barcode(source[required], barcodes[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `barcodes`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `split_by_barcode` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import split_by_barcode

result = split_by_barcode(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    barcodes={"sample1": "AC"}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0226",
  "operation": "split_by_barcode",
  "result_type": "dict",
  "result_preview": "{'sample1': [FastqRecord(identifier='a', sequence='ACGT', quality='IIII')], 'unmatched': [FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 51) split_by_chromosome

- **Purpose:** Document a reusable split by chromosome capability for FASTA and FASTQ workflows. Candidate CAP-0227 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split by chromosome operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.split_by_chromosome(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `split_by_chromosome` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import split_by_chromosome

result = split_by_chromosome(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0227",
  "operation": "split_by_chromosome",
  "result_type": "dict",
  "result_preview": "{'a': [('a', 'ACGT')], 'b': [('b', 'GGAA')]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 52) split_by_length

- **Purpose:** Document a reusable split by length capability for FASTA and FASTQ workflows. Candidate CAP-0228 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split by length operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.split_by_length(source[required], threshold[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `threshold`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `split_by_length` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import split_by_length

result = split_by_length(
    source='>a\nACGT\n>b\nGGAA\n',
    threshold=4
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0228",
  "operation": "split_by_length",
  "result_type": "dict",
  "result_preview": "{'short': [], 'long': [('a', 'ACGT'), ('b', 'GGAA')]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 53) split_by_record

- **Purpose:** Document a reusable split by record capability for FASTA and FASTQ workflows. Candidate CAP-0229 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split by record operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.split_by_record(source[required], records_per_part[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `records_per_part`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `split_by_record` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import split_by_record

result = split_by_record(
    source='>a\nACGT\n>b\nGGAA\n',
    records_per_part=1
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0229",
  "operation": "split_by_record",
  "result_type": "list",
  "result_preview": "[[('a', 'ACGT')], [('b', 'GGAA')]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 54) split_by_size

- **Purpose:** Document a reusable split by size capability for FASTA and FASTQ workflows. Candidate CAP-0230 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split by size operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.split_by_size(source[required], bases_per_part[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `bases_per_part`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `split_by_size` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import split_by_size

result = split_by_size(
    source='>a\nACGT\n>b\nGGAA\n',
    bases_per_part=4
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0230",
  "operation": "split_by_size",
  "result_type": "list",
  "result_preview": "[[('a', 'ACGT')], [('b', 'GGAA')]]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 55) split_by_window

- **Purpose:** Document a reusable split by window capability for FASTA and FASTQ workflows. Candidate CAP-0231 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split by window operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.split_by_window(source[required], window[required], step[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `window`; optional arguments: `step`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `split_by_window` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import split_by_window

result = split_by_window(
    source='>a\nACGT\n>b\nGGAA\n',
    window=3,
    step=2
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0231",
  "operation": "split_by_window",
  "result_type": "list",
  "result_preview": "[('a:0-3', 'ACG'), ('a:2-4', 'GT'), ('b:0-3', 'GGA'), ('b:2-4', 'AA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 56) stats

- **Purpose:** Document a reusable stats capability for FASTA and FASTQ workflows. Candidate CAP-0232 comes from 2 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the stats operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta, fastq). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.stats(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `stats` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import stats

result = stats(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0232",
  "operation": "stats",
  "result_type": "dict",
  "result_preview": "{'records': 2, 'bases': 8, 'lengths': {'count': 2, 'min': 4, 'max': 4, 'mean': 4.0, 'counts': {4: 2}}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 57) subsample

- **Purpose:** Document a reusable subsample capability for FASTA and FASTQ workflows. Candidate CAP-0233 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the subsample operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. Existing-name evidence: `sample`.

- **API:** `dnakit.fasta_fastq.subsample(source[required], fraction[required], seed[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `fraction`; optional arguments: `seed`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `subsample` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import subsample

result = subsample(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    fraction=0.5,
    seed=7
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0233",
  "operation": "subsample",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 58) subsample_records

- **Purpose:** Document a reusable subsample records capability for FASTA and FASTQ workflows. Candidate CAP-0234 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the subsample records operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.subsample_records(source[required], fraction[required], seed[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `fraction`; optional arguments: `seed`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `subsample_records` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import subsample_records

result = subsample_records(
    source='>a\nACGT\n>b\nGGAA\n',
    fraction=0.5,
    seed=7
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0234",
  "operation": "subsample_records",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 59) synchronize_pairs

- **Purpose:** Document a reusable synchronize pairs capability for FASTA and FASTQ workflows. Candidate CAP-0235 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the synchronize pairs operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.synchronize_pairs(first[required], second[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `first`, `second`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `synchronize_pairs` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import synchronize_pairs

result = synchronize_pairs(
    first='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    second='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0235",
  "operation": "synchronize_pairs",
  "result_type": "tuple",
  "result_preview": "([FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')], [FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')])",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 60) trim_adapters

- **Purpose:** Document a reusable trim adapters capability for FASTA and FASTQ workflows. Candidate CAP-0236 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the trim adapters operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.trim_adapters(source[required], adapters[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `adapters`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `trim_adapters` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import trim_adapters

result = trim_adapters(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    adapters=["AAA"]
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0236",
  "operation": "trim_adapters",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 61) trim_low_quality_ends

- **Purpose:** Document a reusable trim low quality ends capability for FASTA and FASTQ workflows. Candidate CAP-0237 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the trim low quality ends operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.trim_low_quality_ends(source[required], minimum_quality[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: `minimum_quality`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `trim_low_quality_ends` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import trim_low_quality_ends

result = trim_low_quality_ends(
    source='@a\nACGT\n+\nIIII\n@b\nGGAA\n+\nIIII\n',
    minimum_quality=20
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0237",
  "operation": "trim_low_quality_ends",
  "result_type": "list",
  "result_preview": "[FastqRecord(identifier='a', sequence='ACGT', quality='IIII'), FastqRecord(identifier='b', sequence='GGAA', quality='IIII')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 62) unmask_regions

- **Purpose:** Document a reusable unmask regions capability for FASTA and FASTQ workflows. Candidate CAP-0238 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the unmask regions operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.unmask_regions(source[required], regions[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`, `regions`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `unmask_regions` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import unmask_regions

result = unmask_regions(
    source='>a\nACGT\n>b\nGGAA\n',
    regions={"a": [(0, 2)]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0238",
  "operation": "unmask_regions",
  "result_type": "list",
  "result_preview": "[('a', 'ACGT'), ('b', 'GGAA')]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 63) unwrap_fasta_lines

- **Purpose:** Document a reusable unwrap fasta lines capability for FASTA and FASTQ workflows. Candidate CAP-0239 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the unwrap fasta lines operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.unwrap_fasta_lines(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `unwrap_fasta_lines` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import unwrap_fasta_lines

result = unwrap_fasta_lines(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0239",
  "operation": "unwrap_fasta_lines",
  "result_type": "str",
  "result_preview": "'>a\\nACGT\\n>b\\nGGAA\\n'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 64) validate

- **Purpose:** Document a reusable validate capability for FASTA and FASTQ workflows. Candidate CAP-0240 comes from 2 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta, fastq). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `generic_operation_name`. Existing-name evidence: `dnakit.standardize.validate|validate|validate_set`.

- **API:** `dnakit.fasta_fastq.validate(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `validate` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import validate

result = validate(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0240",
  "operation": "validate",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'format': 'fasta', 'records': 2, 'duplicate_ids': []}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 65) validate_sequence_ids

- **Purpose:** Document a reusable validate sequence ids capability for FASTA and FASTQ workflows. Candidate CAP-0241 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate sequence ids operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The dependency-free implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. It operates on explicit FASTA/FASTQ records or text and does not invoke an aligner, reference genome, or external database. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.validate_sequence_ids(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `validate_sequence_ids` function returns its native Python value (records, a mapping, or a list of summaries) from explicit FASTA/FASTQ input. It does not infer annotations or call an external backend.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import validate_sequence_ids

result = validate_sequence_ids(
    source='>a\nACGT\n>b\nGGAA\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0241",
  "operation": "validate_sequence_ids",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'ids': ['a', 'b'], 'duplicates': []}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 66) wrap_fasta_lines

- **Purpose:** Document a reusable wrap fasta lines capability for FASTA and FASTQ workflows. Candidate CAP-0242 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the wrap fasta lines operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (fasta). The implementation is available at `src/dnakit/fasta_fastq/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.fasta_fastq.wrap_fasta_lines(source[required], width[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: `width`. FASTA/FASTQ arguments accept text, paths, or materialized records.

- **Output:** The local `wrap_fasta_lines` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.fasta_fastq import wrap_fasta_lines

result = wrap_fasta_lines(
    source='>a\nACGT\n>b\nGGAA\n',
    width=4
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0242",
  "operation": "wrap_fasta_lines",
  "result_type": "str",
  "result_preview": "'>a\\nACGT\\n>b\\nGGAA\\n'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```
