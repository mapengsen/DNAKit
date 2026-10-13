# Alignment

This category contains **34** documented candidate capabilities; 34 have a locally callable formal entry (5 partial due to binary-format dependencies; 0 bounded external backends; 1 dependency-adapter only). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) align_long_reads

- **Purpose:** Document a reusable align long reads capability for Alignment workflows. Candidate CAP-0001 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the align long reads operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.align_long_reads(input[required], backend[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `backend`. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `align_long_reads` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import align_long_reads

result = align_long_reads(
    input={"reads": ["ACGT"], "reference": "ACGTACGT"}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0001",
  "operation": "align_long_reads",
  "result_type": "list",
  "result_preview": "[{'query_name': 'read_1', 'flag': 0, 'reference': 'reference', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'backend': 'naive'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 2) align_paired_reads

- **Purpose:** Document a reusable align paired reads capability for Alignment workflows. Candidate CAP-0002 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the align paired reads operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.align_paired_reads(input[required], backend[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `backend`. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `align_paired_reads` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import align_paired_reads

result = align_paired_reads(
    input={"reads": ["ACGT"], "reference": "ACGTACGT"}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0002",
  "operation": "align_paired_reads",
  "result_type": "list",
  "result_preview": "[{'query_name': 'read_1', 'flag': 0, 'reference': 'reference', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'backend': 'naive'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 3) align_reads

- **Purpose:** Document a reusable align reads capability for Alignment workflows. Candidate CAP-0003 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the align reads operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.align_reads(input[required], backend[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `backend`. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `align_reads` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import align_reads

result = align_reads(
    input={"reads": ["ACGT"], "reference": "ACGTACGT"}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0003",
  "operation": "align_reads",
  "result_type": "list",
  "result_preview": "[{'query_name': 'read_1', 'flag': 0, 'reference': 'reference', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'backend': 'naive'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 4) align_to_graph

- **Purpose:** Document a reusable align to graph capability for Alignment workflows. Candidate CAP-0004 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the align to graph operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.align_to_graph(input[required], backend[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `backend`. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `align_to_graph` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import align_to_graph

result = align_to_graph(
    input={"reads": ["ACGT"], "reference": "ACGTACGT"}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0004",
  "operation": "align_to_graph",
  "result_type": "list",
  "result_preview": "[{'query_name': 'read_1', 'flag': 0, 'reference': 'reference', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'backend': 'naive'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 5) bam_to_cram

- **Purpose:** Document a reusable bam to cram capability for Alignment workflows. Candidate CAP-0005 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the bam to cram operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/alignment/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This entry is partial: the local example handles SAM text; binary BAM/CRAM requires a validated pysam or equivalent backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.bam_to_cram(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `bam_to_cram` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview. Binary BAM/CRAM parsing or conversion is not validated and remains pending.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import bam_to_cram

result = bam_to_cram(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_partial",
  "capability_id": "CAP-0005",
  "operation": "bam_to_cram",
  "result_type": "bytes",
  "result_preview": "<bytes length=34>",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_partial"
  }
}
```

## 6) bam_to_sam

- **Purpose:** Document a reusable bam to sam capability for Alignment workflows. Candidate CAP-0006 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the bam to sam operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/alignment/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This entry is partial: the local example handles SAM text; binary BAM/CRAM requires a validated pysam or equivalent backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.bam_to_sam(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `bam_to_sam` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview. Binary BAM/CRAM parsing or conversion is not validated and remains pending.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import bam_to_sam

result = bam_to_sam(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_partial",
  "capability_id": "CAP-0006",
  "operation": "bam_to_sam",
  "result_type": "str",
  "result_preview": "'r1\\t0\\tchr1\\t1\\t60\\t4M\\t*\\t0\\t0\\tACGT\\tIIII\\n'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_partial"
  }
}
```

## 7) build_index

- **Purpose:** Document a reusable build index capability for Alignment workflows. Candidate CAP-0007 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the build index operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `generic_operation_name`. Existing-name evidence: `build_index|dnakit.download.build_index`.

- **API:** `dnakit.alignment.build_index(input[required], k[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `k`. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `build_index` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import build_index

result = build_index(
    input={"reads": ["ACGT"], "reference": "ACGTACGT"},
    k=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0007",
  "operation": "build_index",
  "result_type": "dict",
  "result_preview": "{'ACG': [0, 4], 'CGT': [1, 5], 'GTA': [2], 'TAC': [3]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 8) calculate_base_depth

- **Purpose:** Document a reusable calculate base depth capability for Alignment workflows. Candidate CAP-0008 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate base depth operation, and return structured results with provenance and diagnostics. The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.calculate_base_depth(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `calculate_base_depth` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import calculate_base_depth

result = calculate_base_depth(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0008",
  "operation": "calculate_base_depth",
  "result_type": "dict",
  "result_preview": "{'chr1': {0: 1, 1: 1, 2: 1, 3: 1}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 9) calculate_clipping_profile

- **Purpose:** Document a reusable calculate clipping profile capability for Alignment workflows. Candidate CAP-0009 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate clipping profile operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/alignment/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.calculate_clipping_profile(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local adapter `calculate_clipping_profile` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import calculate_clipping_profile

result = calculate_clipping_profile(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0009",
  "operation": "calculate_clipping_profile",
  "result_type": "dict",
  "result_preview": "{'module': 'alignment', 'operation': 'calculate_clipping_profile', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 45, 'sequence_bases': 14, 'lengths': [14], 'sha256': 'b20b0e08c38e755055c23270535b5a3daa4dda96fa1...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 10) calculate_coverage

- **Purpose:** Document a reusable calculate coverage capability for Alignment workflows. Candidate CAP-0010 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate coverage operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.calculate_coverage(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `calculate_coverage` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import calculate_coverage

result = calculate_coverage(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0010",
  "operation": "calculate_coverage",
  "result_type": "dict",
  "result_preview": "{'covered_bases': 4, 'total_depth': 4, 'mean_depth': 1.0, 'by_reference': {'chr1': {0: 1, 1: 1, 2: 1, 3: 1}}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 11) calculate_depth

- **Purpose:** Document a reusable calculate depth capability for Alignment workflows. Candidate CAP-0011 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate depth operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.calculate_depth(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `calculate_depth` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import calculate_depth

result = calculate_depth(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0011",
  "operation": "calculate_depth",
  "result_type": "dict",
  "result_preview": "{'covered_bases': 4, 'total_depth': 4, 'mean_depth': 1.0, 'by_reference': {'chr1': {0: 1, 1: 1, 2: 1, 3: 1}}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 12) calculate_insert_size

- **Purpose:** Document a reusable calculate insert size capability for Alignment workflows. Candidate CAP-0012 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate insert size operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.calculate_insert_size(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `calculate_insert_size` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import calculate_insert_size

result = calculate_insert_size(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0012",
  "operation": "calculate_insert_size",
  "result_type": "dict",
  "result_preview": "{'count': 0, 'mean': 0.0, 'values': []}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 13) calculate_mapping_rate

- **Purpose:** Document a reusable calculate mapping rate capability for Alignment workflows. Candidate CAP-0013 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate mapping rate operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.calculate_mapping_rate(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `calculate_mapping_rate` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import calculate_mapping_rate

result = calculate_mapping_rate(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0013",
  "operation": "calculate_mapping_rate",
  "result_type": "float",
  "result_preview": "1.0",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 14) calculate_mismatch_rate

- **Purpose:** Document a reusable calculate mismatch rate capability for Alignment workflows. Candidate CAP-0014 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate mismatch rate operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.calculate_mismatch_rate(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `calculate_mismatch_rate` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import calculate_mismatch_rate

result = calculate_mismatch_rate(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0014",
  "operation": "calculate_mismatch_rate",
  "result_type": "float",
  "result_preview": "0.0",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 15) calculate_read_group_stats

- **Purpose:** Document a reusable calculate read group stats capability for Alignment workflows. Candidate CAP-0015 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate read group stats operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.calculate_read_group_stats(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `calculate_read_group_stats` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import calculate_read_group_stats

result = calculate_read_group_stats(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0015",
  "operation": "calculate_read_group_stats",
  "result_type": "dict",
  "result_preview": "{'unknown': {'records': 1, 'mapped': 1}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 16) calculate_soft_clipping

- **Purpose:** Document a reusable calculate soft clipping capability for Alignment workflows. Candidate CAP-0016 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate soft clipping operation, and return structured results with provenance and diagnostics. The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.calculate_soft_clipping(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `calculate_soft_clipping` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import calculate_soft_clipping

result = calculate_soft_clipping(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0016",
  "operation": "calculate_soft_clipping",
  "result_type": "dict",
  "result_preview": "{'r1': 0}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 17) calculate_strand_balance

- **Purpose:** Document a reusable calculate strand balance capability for Alignment workflows. Candidate CAP-0017 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate strand balance operation, and return structured results with provenance and diagnostics. The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.calculate_strand_balance(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `calculate_strand_balance` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import calculate_strand_balance

result = calculate_strand_balance(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0017",
  "operation": "calculate_strand_balance",
  "result_type": "dict",
  "result_preview": "{'forward': 1, 'reverse': 0, 'forward_fraction': 1.0}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 18) calculate_strand_bias

- **Purpose:** Document a reusable calculate strand bias capability for Alignment workflows. Candidate CAP-0018 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate strand bias operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.calculate_strand_bias(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `calculate_strand_bias` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import calculate_strand_bias

result = calculate_strand_bias(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0018",
  "operation": "calculate_strand_bias",
  "result_type": "dict",
  "result_preview": "{'bias': 1.0, 'forward': 1, 'reverse': 0, 'forward_fraction': 1.0}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 19) call_consensus_from_pileup

- **Purpose:** Document a reusable call consensus from pileup capability for Alignment workflows. Candidate CAP-0019 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the call consensus from pileup operation, and return structured results with provenance and diagnostics. The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.call_consensus_from_pileup(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `call_consensus_from_pileup` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import call_consensus_from_pileup

result = call_consensus_from_pileup(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0019",
  "operation": "call_consensus_from_pileup",
  "result_type": "dict",
  "result_preview": "{'chr1': 'ACGT'}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 20) filter_bam_by_flags

- **Purpose:** Document a reusable filter bam by flags capability for Alignment workflows. Candidate CAP-0020 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter bam by flags operation, and return structured results with provenance and diagnostics. The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.filter_bam_by_flags(input[required], include_flags[optional], exclude_flags[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `include_flags`, `exclude_flags`. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `filter_bam_by_flags` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import filter_bam_by_flags

result = filter_bam_by_flags(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0020",
  "operation": "filter_bam_by_flags",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 21) filter_bam_by_mapping_quality

- **Purpose:** Document a reusable filter bam by mapping quality capability for Alignment workflows. Candidate CAP-0021 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter bam by mapping quality operation, and return structured results with provenance and diagnostics. The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.filter_bam_by_mapping_quality(input[required], minimum[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `minimum`. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `filter_bam_by_mapping_quality` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import filter_bam_by_mapping_quality

result = filter_bam_by_mapping_quality(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0021",
  "operation": "filter_bam_by_mapping_quality",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 22) generate_pileup

- **Purpose:** Document a reusable generate pileup capability for Alignment workflows. Candidate CAP-0022 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the generate pileup operation, and return structured results with provenance and diagnostics. The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.generate_pileup(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `generate_pileup` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import generate_pileup

result = generate_pileup(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0022",
  "operation": "generate_pileup",
  "result_type": "list",
  "result_preview": "[{'reference': 'chr1', 'position': 0, 'depth': 1}, {'reference': 'chr1', 'position': 1, 'depth': 1}, {'reference': 'chr1', 'position': 2, 'depth': 1}, {'reference': 'chr1', 'position': 3, 'depth': 1}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 23) index_alignments

- **Purpose:** Document a reusable index alignments capability for Alignment workflows. Candidate CAP-0023 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the index alignments operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.index_alignments(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `index_alignments` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import index_alignments

result = index_alignments(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0023",
  "operation": "index_alignments",
  "result_type": "dict",
  "result_preview": "{'chr1': [0]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 24) mark_duplicates

- **Purpose:** Document a reusable mark duplicates capability for Alignment workflows. Candidate CAP-0024 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the mark duplicates operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.mark_duplicates(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `mark_duplicates` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import mark_duplicates

result = mark_duplicates(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0024",
  "operation": "mark_duplicates",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 25) merge_bam_files

- **Purpose:** Document a reusable merge bam files capability for Alignment workflows. Candidate CAP-0025 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge bam files operation, and return structured results with provenance and diagnostics. The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.merge_bam_files(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `merge_bam_files` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import merge_bam_files

result = merge_bam_files(
    input=['@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n']
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0025",
  "operation": "merge_bam_files",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 26) parse_bam_records

- **Purpose:** Document a reusable parse bam records capability for Alignment workflows. Candidate CAP-0026 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the parse bam records operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/alignment/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This entry is partial: the local example handles SAM text; binary BAM/CRAM requires a validated pysam or equivalent backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.parse_bam_records(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `parse_bam_records` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview. Binary BAM/CRAM parsing or conversion is not validated and remains pending.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import parse_bam_records

result = parse_bam_records(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_partial",
  "capability_id": "CAP-0026",
  "operation": "parse_bam_records",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_partial"
  }
}
```

## 27) parse_cram_records

- **Purpose:** Document a reusable parse cram records capability for Alignment workflows. Candidate CAP-0027 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the parse cram records operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/alignment/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This entry is partial: the local example handles SAM text; binary BAM/CRAM requires a validated pysam or equivalent backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.parse_cram_records(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `parse_cram_records` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview. Binary BAM/CRAM parsing or conversion is not validated and remains pending.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import parse_cram_records

result = parse_cram_records(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_partial",
  "capability_id": "CAP-0027",
  "operation": "parse_cram_records",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_partial"
  }
}
```

## 28) parse_sam_records

- **Purpose:** Document a reusable parse sam records capability for Alignment workflows. Candidate CAP-0028 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the parse sam records operation, and return structured results with provenance and diagnostics. The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.parse_sam_records(source[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `source`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `parse_sam_records` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import parse_sam_records

result = parse_sam_records(
    source='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0028",
  "operation": "parse_sam_records",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 29) plot_alignment_view

- **Purpose:** Document a reusable plot alignment view capability for Alignment workflows. Candidate CAP-0029 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the plot alignment view operation, and return structured results with provenance and diagnostics. The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `no_exact_public_match`. Existing-name evidence: `plot_alignment`.

- **API:** `dnakit.alignment.plot_alignment_view(input[required], width[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `width`. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `plot_alignment_view` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import plot_alignment_view

result = plot_alignment_view(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n',
    width=3
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0029",
  "operation": "plot_alignment_view",
  "result_type": "str",
  "result_preview": "'<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"24\" height=\"30\"><text x=\"0\" y=\"18\">r1 chr1:1 ACGT</text></svg>'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 30) plot_coverage_track

- **Purpose:** Document a reusable plot coverage track capability for Alignment workflows. Candidate CAP-0030 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the plot coverage track operation, and return structured results with provenance and diagnostics. The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.plot_coverage_track(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `plot_coverage_track` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import plot_coverage_track

result = plot_coverage_track(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0030",
  "operation": "plot_coverage_track",
  "result_type": "str",
  "result_preview": "'<svg xmlns=\"http://www.w3.org/2000/svg\"><text>mean_depth=1.0</text></svg>'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 31) remove_duplicates

- **Purpose:** Document a reusable remove duplicates capability for Alignment workflows. Candidate CAP-0031 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the remove duplicates operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.remove_duplicates(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `remove_duplicates` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import remove_duplicates

result = remove_duplicates(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0031",
  "operation": "remove_duplicates",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 32) sam_to_bam

- **Purpose:** Document a reusable sam to bam capability for Alignment workflows. Candidate CAP-0032 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sam to bam operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/alignment/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This entry is partial: the local example handles SAM text; binary BAM/CRAM requires a validated pysam or equivalent backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.sam_to_bam(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `sam_to_bam` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview. Binary BAM/CRAM parsing or conversion is not validated and remains pending.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import sam_to_bam

result = sam_to_bam(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_partial",
  "capability_id": "CAP-0032",
  "operation": "sam_to_bam",
  "result_type": "bytes",
  "result_preview": "<bytes length=34>",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_partial"
  }
}
```

## 33) sort_alignments

- **Purpose:** Document a reusable sort alignments capability for Alignment workflows. Candidate CAP-0033 comes from 10 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sort alignments operation, and return structured results with provenance and diagnostics. The optional `backend` parameter selects one supported variant (bowtie2, bwa, graphmap, hisat2, minimap2, ngmlr, pbmm2, star, vg, winnowmap). The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `merged_dimension`. Existing-name evidence: `plot_alignment`.

- **API:** `dnakit.alignment.sort_alignments(input[required], by[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `by`. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `sort_alignments` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import sort_alignments

result = sort_alignments(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0033",
  "operation": "sort_alignments",
  "result_type": "list",
  "result_preview": "[{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 34) split_bam_by_region

- **Purpose:** Document a reusable split bam by region capability for Alignment workflows. Candidate CAP-0034 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split bam by region operation, and return structured results with provenance and diagnostics. The deterministic local implementation is available at `src/dnakit/alignment/_candidates.py` and is exported through the formal DNAKit namespace. It operates on SAM-like text or materialized records; binary BAM/CRAM parsing, reference-indexed alignment, and external aligners remain outside this boundary. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.alignment.split_bam_by_region(input[required], regions[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`, `regions`; optional arguments: no. Alignment input accepts SAM text, record mappings, or the local read/reference envelope.

- **Output:** The local `split_bam_by_region` function returns its deterministic Python result from SAM-like text or materialized records. It does not parse binary BAM/CRAM or invoke a reference-indexed external aligner.

- **Example code:**

```python
from pathlib import Path

from dnakit.alignment import split_bam_by_region

result = split_bam_by_region(
    input='@HD\tVN:1.6\nr1\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n',
    regions={"chr1": [(0, 4)]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0034",
  "operation": "split_bam_by_region",
  "result_type": "dict",
  "result_preview": "{'chr1': [{'query_name': 'r1', 'flag': 0, 'reference': 'chr1', 'position': 0, 'mapping_quality': 60, 'cigar': '4M', 'sequence': 'ACGT', 'quality': 'IIII', 'tags': {}}]}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```
