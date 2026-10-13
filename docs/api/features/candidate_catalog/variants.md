# Variant analysis

This category contains **63** documented candidate capabilities; 63 have a locally callable formal entry (0 partial due to binary-format dependencies; 0 bounded external backends; 48 dependency-adapter only). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) annotate

- **Purpose:** Document a reusable annotate capability for Variant analysis workflows. Candidate CAP-0759 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the annotate operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.annotate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `annotate` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import annotate

result = annotate(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0759",
  "operation": "annotate",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'annotate', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a21887...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 2) annotate_coding_effect

- **Purpose:** Document a reusable annotate coding effect capability for Variant analysis workflows. Candidate CAP-0760 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the annotate coding effect operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.annotate_coding_effect(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `annotate_coding_effect` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import annotate_coding_effect

result = annotate_coding_effect(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0760",
  "operation": "annotate_coding_effect",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'annotate_coding_effect', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee67...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 3) annotate_conservation_score

- **Purpose:** Document a reusable annotate conservation score capability for Variant analysis workflows. Candidate CAP-0761 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the annotate conservation score operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.annotate_conservation_score(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `annotate_conservation_score` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import annotate_conservation_score

result = annotate_conservation_score(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0761",
  "operation": "annotate_conservation_score",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'annotate_conservation_score', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d9...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 4) annotate_mt

- **Purpose:** Document a reusable annotate mt capability for Variant analysis workflows. Candidate CAP-0762 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the annotate mt operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.annotate_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `annotate_mt` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import annotate_mt

result = annotate_mt(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0762",
  "operation": "annotate_mt",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'annotate_mt', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a21...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 5) annotate_pathogenicity_score

- **Purpose:** Document a reusable annotate pathogenicity score capability for Variant analysis workflows. Candidate CAP-0763 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the annotate pathogenicity score operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.annotate_pathogenicity_score(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `annotate_pathogenicity_score` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import annotate_pathogenicity_score

result = annotate_pathogenicity_score(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0763",
  "operation": "annotate_pathogenicity_score",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'annotate_pathogenicity_score', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 6) annotate_population_frequency

- **Purpose:** Document a reusable annotate population frequency capability for Variant analysis workflows. Candidate CAP-0764 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the annotate population frequency operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.annotate_population_frequency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `annotate_population_frequency` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import annotate_population_frequency

result = annotate_population_frequency(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0764",
  "operation": "annotate_population_frequency",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'annotate_population_frequency', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 7) annotate_regulatory_effect

- **Purpose:** Document a reusable annotate regulatory effect capability for Variant analysis workflows. Candidate CAP-0765 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the annotate regulatory effect operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.annotate_regulatory_effect(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `annotate_regulatory_effect` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import annotate_regulatory_effect

result = annotate_regulatory_effect(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0765",
  "operation": "annotate_regulatory_effect",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'annotate_regulatory_effect', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 8) annotate_splice_effect

- **Purpose:** Document a reusable annotate splice effect capability for Variant analysis workflows. Candidate CAP-0766 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the annotate splice effect operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.annotate_splice_effect(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `annotate_splice_effect` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import annotate_splice_effect

result = annotate_splice_effect(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0766",
  "operation": "annotate_splice_effect",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'annotate_splice_effect', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee67...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 9) annotate_variant_consequence

- **Purpose:** Document a reusable annotate variant consequence capability for Variant analysis workflows. Candidate CAP-0767 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the annotate variant consequence operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.annotate_variant_consequence(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `annotate_variant_consequence` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import annotate_variant_consequence

result = annotate_variant_consequence(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0767",
  "operation": "annotate_variant_consequence",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'annotate_variant_consequence', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 10) calculate_allele_frequency

- **Purpose:** Document a reusable calculate allele frequency capability for Variant analysis workflows. Candidate CAP-0768 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate allele frequency operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_allele_frequency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics.

- **Output:** The local `calculate_allele_frequency` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_allele_frequency

result = calculate_allele_frequency(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0768",
  "operation": "calculate_allele_frequency",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_allele_frequency', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': ['G'], 'alt_count': 3, 'called_alleles': 6, 'frequency': 0.5}], 'total_alt_alleles': 3, 'total_called_alleles': 6, 'aggregate_frequency': 0.5}, 'provenance': ...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 11) calculate_copy_number_profile

- **Purpose:** Document a reusable calculate copy number profile capability for Variant analysis workflows. Candidate CAP-0769 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate copy number profile operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_copy_number_profile(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_copy_number_profile` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_copy_number_profile

result = calculate_copy_number_profile(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0769",
  "operation": "calculate_copy_number_profile",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'calculate_copy_number_profile', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 12) calculate_genotype_likelihood

- **Purpose:** Document a reusable calculate genotype likelihood capability for Variant analysis workflows. Candidate CAP-0770 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate genotype likelihood operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_genotype_likelihood(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_genotype_likelihood` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_genotype_likelihood

result = calculate_genotype_likelihood(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0770",
  "operation": "calculate_genotype_likelihood",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'calculate_genotype_likelihood', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 13) calculate_haplotype_frequency

- **Purpose:** Document a reusable calculate haplotype frequency capability for Variant analysis workflows. Candidate CAP-0771 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate haplotype frequency operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_haplotype_frequency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_haplotype_frequency` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_haplotype_frequency

result = calculate_haplotype_frequency(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0771",
  "operation": "calculate_haplotype_frequency",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'calculate_haplotype_frequency', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 14) calculate_heterozygosity

- **Purpose:** Document a reusable calculate heterozygosity capability for Variant analysis workflows. Candidate CAP-0772 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate heterozygosity operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_heterozygosity(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics.

- **Output:** The local `calculate_heterozygosity` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_heterozygosity

result = calculate_heterozygosity(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0772",
  "operation": "calculate_heterozygosity",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_heterozygosity', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': ['G'], 'heterozygous_samples': 1, 'called_samples': 3, 'observed_heterozygosity': 0.3333333333333333}], 'heterozygous_samples': 1, 'called_samples': 3, 'observe...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 15) calculate_hwe_pvalue

- **Purpose:** Document a reusable calculate hwe pvalue capability for Variant analysis workflows. Candidate CAP-0773 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate hwe pvalue operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_hwe_pvalue(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_hwe_pvalue` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_hwe_pvalue

result = calculate_hwe_pvalue(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0773",
  "operation": "calculate_hwe_pvalue",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'calculate_hwe_pvalue', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 16) calculate_identity_by_descent

- **Purpose:** Document a reusable calculate identity by descent capability for Variant analysis workflows. Candidate CAP-0774 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate identity by descent operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_identity_by_descent(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_identity_by_descent` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_identity_by_descent

result = calculate_identity_by_descent(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0774",
  "operation": "calculate_identity_by_descent",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'calculate_identity_by_descent', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 17) calculate_ld_matrix

- **Purpose:** Document a reusable calculate ld matrix capability for Variant analysis workflows. Candidate CAP-0775 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate ld matrix operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_ld_matrix(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_ld_matrix` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_ld_matrix

result = calculate_ld_matrix(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0775",
  "operation": "calculate_ld_matrix",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'calculate_ld_matrix', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e0...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 18) calculate_minor_allele_frequency

- **Purpose:** Document a reusable calculate minor allele frequency capability for Variant analysis workflows. Candidate CAP-0776 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate minor allele frequency operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_minor_allele_frequency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics. Explicit `genotypes` are required; no missing-genotype imputation is performed.

- **Output:** The local `calculate_minor_allele_frequency` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_minor_allele_frequency

result = calculate_minor_allele_frequency(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0776",
  "operation": "calculate_minor_allele_frequency",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_minor_allele_frequency', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': ['G'], 'called_alleles': 6, 'allele_frequencies': {'0': 0.5, '1': 0.5}, 'minor_allele_index': 0, 'minor_allele': 'A', 'maf': 0.5}], 'variant_count': 1, ...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 19) calculate_missingness

- **Purpose:** Document a reusable calculate missingness capability for Variant analysis workflows. Candidate CAP-0777 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate missingness operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_missingness(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics.

- **Output:** The local `calculate_missingness` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_missingness

result = calculate_missingness(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0777",
  "operation": "calculate_missingness",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_missingness', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': ['G'], 'missing_samples': 0, 'sample_count': 3, 'missingness': 0.0}], 'missing_samples': 0, 'sample_count': 3, 'missingness': 0.0}, 'provenance': {'implementation'...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 20) calculate_runs_of_homozygosity

- **Purpose:** Document a reusable calculate runs of homozygosity capability for Variant analysis workflows. Candidate CAP-0778 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate runs of homozygosity operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_runs_of_homozygosity(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `calculate_runs_of_homozygosity` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_runs_of_homozygosity

result = calculate_runs_of_homozygosity(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0778",
  "operation": "calculate_runs_of_homozygosity",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'calculate_runs_of_homozygosity', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf42959...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 21) calculate_transition_transversion

- **Purpose:** Document a reusable calculate transition transversion capability for Variant analysis workflows. Candidate CAP-0779 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate transition transversion operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_transition_transversion(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics.

- **Output:** The local `calculate_transition_transversion` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_transition_transversion

result = calculate_transition_transversion(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0779",
  "operation": "calculate_transition_transversion",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_transition_transversion', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'classification': 'transition'}], 'transitions': 1, 'transversions': 0, 'excluded': 0, 'ti_tv_ratio': None}, 'provenance': {'implementation': 'dna...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 22) calculate_variant_density

- **Purpose:** Document a reusable calculate variant density capability for Variant analysis workflows. Candidate CAP-0780 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the calculate variant density operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.calculate_variant_density(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics. `options.region_length` is required for variant density.

- **Output:** The local `calculate_variant_density` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import calculate_variant_density

result = calculate_variant_density(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]},
    options={"region_length": 100}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0780",
  "operation": "calculate_variant_density",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'calculate_variant_density', 'result': {'variant_count': 1, 'region_length': 100, 'density_per_base': 0.01, 'density_per_megabase': 10000.0, 'by_chromosome': {'chr1': 1}}, 'provenance': {'implementation': 'dnakit.native_variant_statistics', 'backend': 'local_d...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 23) call

- **Purpose:** Document a reusable call capability for Variant analysis workflows. Candidate CAP-0781 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the call operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `generic_operation_name`. Existing-name evidence: `call|dnakit.databases.call`.

- **API:** `dnakit.variants.call(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `call` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import call

result = call(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0781",
  "operation": "call",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'call', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a21887'}, ...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 24) call_mt

- **Purpose:** Document a reusable call mt capability for Variant analysis workflows. Candidate CAP-0782 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the call mt operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.call_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `call_mt` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import call_mt

result = call_mt(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0782",
  "operation": "call_mt",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'call_mt', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a21887'...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 25) compare

- **Purpose:** Document a reusable compare capability for Variant analysis workflows. Candidate CAP-0783 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compare operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `generic_operation_name`. Existing-name evidence: `compare|dnakit.similarity.compare`.

- **API:** `dnakit.variants.compare(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping with `left` and `right` variant collections; matching uses exact CHROM/POS/REF/ALT identity without normalization.

- **Output:** The local `compare` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import compare

result = compare(
    input={"left": {"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"}]}, "right": {"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G"}, {"chrom": "chr1", "pos": 20, "ref": "C", "alt": "T"}]}}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0783",
  "operation": "compare",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'compare', 'result': {'left_record_count': 1, 'right_record_count': 2, 'left_allele_count': 1, 'right_allele_count': 2, 'shared_allele_count': 1, 'only_left': [], 'only_right': [{'chrom': 'chr1', 'pos': 20, 'ref': 'C', 'alt': 'T'}], 'shared': [{'chrom': 'chr1'...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 26) compare_mt

- **Purpose:** Document a reusable compare mt capability for Variant analysis workflows. Candidate CAP-0784 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compare mt operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. Existing-name evidence: `compare`.

- **API:** `dnakit.variants.compare_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compare_mt` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import compare_mt

result = compare_mt(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0784",
  "operation": "compare_mt",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'compare_mt', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a218...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 27) compare_vcfs

- **Purpose:** Document a reusable compare vcfs capability for Variant analysis workflows. Candidate CAP-0785 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compare vcfs operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.compare_vcfs(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compare_vcfs` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import compare_vcfs

result = compare_vcfs(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0785",
  "operation": "compare_vcfs",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'compare_vcfs', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 28) construct_haplotypes

- **Purpose:** Document a reusable construct haplotypes capability for Variant analysis workflows. Candidate CAP-0786 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the construct haplotypes operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.construct_haplotypes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `construct_haplotypes` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import construct_haplotypes

result = construct_haplotypes(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0786",
  "operation": "construct_haplotypes",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'construct_haplotypes', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 29) count

- **Purpose:** Document a reusable count capability for Variant analysis workflows. Candidate CAP-0787 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the count operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.count(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a sequence or mapping containing VCF-like records; counts are based on explicit CHROM/POS/REF/ALT fields.

- **Output:** The local `count` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import count

result = count(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0787",
  "operation": "count",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'count', 'result': {'record_count': 1, 'alternate_allele_count': 1, 'by_chromosome': {'chr1': 1}, 'allele_type_counts': {'snp': 1}}, 'provenance': {'implementation': 'dnakit.native_variant_statistics', 'backend': 'local_dependency_free', 'coordinate_system': '...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 30) count_mt

- **Purpose:** Document a reusable count mt capability for Variant analysis workflows. Candidate CAP-0788 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the count mt operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.count_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `count_mt` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import count_mt

result = count_mt(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0788",
  "operation": "count_mt",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'count_mt', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a21887...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 31) detect_breakpoints

- **Purpose:** Document a reusable detect breakpoints capability for Variant analysis workflows. Candidate CAP-0789 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect breakpoints operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.detect_breakpoints(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `detect_breakpoints` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import detect_breakpoints

result = detect_breakpoints(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0789",
  "operation": "detect_breakpoints",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'detect_breakpoints', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 32) detect_low_frequency_variants

- **Purpose:** Document a reusable detect low frequency variants capability for Variant analysis workflows. Candidate CAP-0790 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect low frequency variants operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.detect_low_frequency_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `detect_low_frequency_variants` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import detect_low_frequency_variants

result = detect_low_frequency_variants(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0790",
  "operation": "detect_low_frequency_variants",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'detect_low_frequency_variants', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 33) detect_mobile_element_insertions

- **Purpose:** Document a reusable detect mobile element insertions capability for Variant analysis workflows. Candidate CAP-0791 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect mobile element insertions operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.detect_mobile_element_insertions(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `detect_mobile_element_insertions` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import detect_mobile_element_insertions

result = detect_mobile_element_insertions(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0791",
  "operation": "detect_mobile_element_insertions",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'detect_mobile_element_insertions', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 34) detect_mosaic_variants

- **Purpose:** Document a reusable detect mosaic variants capability for Variant analysis workflows. Candidate CAP-0792 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect mosaic variants operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.detect_mosaic_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `detect_mosaic_variants` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import detect_mosaic_variants

result = detect_mosaic_variants(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0792",
  "operation": "detect_mosaic_variants",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'detect_mosaic_variants', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee67...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 35) detect_repeat_expansions

- **Purpose:** Document a reusable detect repeat expansions capability for Variant analysis workflows. Candidate CAP-0793 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect repeat expansions operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.detect_repeat_expansions(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `detect_repeat_expansions` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import detect_repeat_expansions

result = detect_repeat_expansions(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0793",
  "operation": "detect_repeat_expansions",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'detect_repeat_expansions', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 36) filter

- **Purpose:** Document a reusable filter capability for Variant analysis workflows. Candidate CAP-0794 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.filter(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `filter` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import filter

result = filter(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0794",
  "operation": "filter",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'filter', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a21887'}...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 37) filter_mt

- **Purpose:** Document a reusable filter mt capability for Variant analysis workflows. Candidate CAP-0795 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter mt operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.filter_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `filter_mt` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import filter_mt

result = filter_mt(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0795",
  "operation": "filter_mt",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'filter_mt', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a2188...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 38) filter_variants_by_allele_frequency

- **Purpose:** Document a reusable filter variants by allele frequency capability for Variant analysis workflows. Candidate CAP-0796 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter variants by allele frequency operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.filter_variants_by_allele_frequency(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref` and `alt`, with explicit `genotypes` for genotype-based statistics. Optional `options.min_frequency`, `max_frequency`, and `allele_index` select the explicit frequency range.

- **Output:** The local `filter_variants_by_allele_frequency` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_allele_frequency

result = filter_variants_by_allele_frequency(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]},
    options={"min_frequency": 0.1, "max_frequency": 0.9}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0796",
  "operation": "filter_variants_by_allele_frequency",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_allele_frequency', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'genotypes': ['0/1', '1/1', '0/0']}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'allele_index': 1, 'frequency': 0.5, 'kept': True}], 'input_recor...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 39) filter_variants_by_depth

- **Purpose:** Document a reusable filter variants by depth capability for Variant analysis workflows. Candidate CAP-0797 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter variants by depth operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.filter_variants_by_depth(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref`, `alt`, and an explicit `DP`/`depth` value. `options.min_depth` is required; optional `max_depth` and `keep_missing` bound explicit DP values.

- **Output:** The local `filter_variants_by_depth` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_depth

result = filter_variants_by_depth(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "DP": 30}]},
    options={"min_depth": 10, "max_depth": 100}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0797",
  "operation": "filter_variants_by_depth",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_depth', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'DP': 30}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'depth': 30.0, 'kept': True}], 'input_record_count': 1, 'output_record_count': 1, 'removed_record_coun...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 40) filter_variants_by_genotype_quality

- **Purpose:** Document a reusable filter variants by genotype quality capability for Variant analysis workflows. Candidate CAP-0798 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter variants by genotype quality operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.filter_variants_by_genotype_quality(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires explicit `GQ`, `genotype_quality`, `genotype_qualities`, or sample GQ values. `options.min_genotype_quality` or `min_gq` is required; GQ values are taken from explicit records or samples.

- **Output:** The local `filter_variants_by_genotype_quality` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_genotype_quality

result = filter_variants_by_genotype_quality(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "GQ": 45}]},
    options={"min_genotype_quality": 20}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0798",
  "operation": "filter_variants_by_genotype_quality",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_genotype_quality', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'GQ': 45}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'genotype_qualities': [45.0], 'minimum_observed': 45.0, 'kept': True}], 'input_record_count...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 41) filter_variants_by_missingness

- **Purpose:** Document a reusable filter variants by missingness capability for Variant analysis workflows. Candidate CAP-0799 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter variants by missingness operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.filter_variants_by_missingness(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref`, `alt`, and either explicit `genotypes` or `missing_count/sample_count`. `options.max_missingness` is required and must be between 0 and 1.

- **Output:** The local `filter_variants_by_missingness` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_missingness

result = filter_variants_by_missingness(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "./.", "0/0"]}]},
    options={"max_missingness": 0.5}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0799",
  "operation": "filter_variants_by_missingness",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_missingness', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'genotypes': ['0/1', './.', '0/0']}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'missingness': 0.3333333333333333, 'kept': True}], 'input_record_count...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 42) filter_variants_by_quality

- **Purpose:** Document a reusable filter variants by quality capability for Variant analysis workflows. Candidate CAP-0800 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter variants by quality operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.filter_variants_by_quality(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a mapping containing VCF-like `variants`; each record requires `chrom`, positive 1-based `pos`, `ref`, `alt`, and a finite `qual`/`QUAL` value. `options.min_quality` is required; missing QUAL values are excluded unless `keep_missing=True`.

- **Output:** The local `filter_variants_by_quality` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import filter_variants_by_quality

result = filter_variants_by_quality(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "qual": 50}]},
    options={"min_quality": 20}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0800",
  "operation": "filter_variants_by_quality",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'filter_variants_by_quality', 'result': {'variants': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'qual': 50}], 'metrics': [{'chrom': 'chr1', 'pos': 10, 'quality': 50.0, 'kept': True}], 'input_record_count': 1, 'output_record_count': 1, 'removed_recor...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 43) genotype_structural_variants

- **Purpose:** Document a reusable genotype structural variants capability for Variant analysis workflows. Candidate CAP-0801 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the genotype structural variants operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.genotype_structural_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `genotype_structural_variants` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import genotype_structural_variants

result = genotype_structural_variants(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0801",
  "operation": "genotype_structural_variants",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'genotype_structural_variants', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 44) impute_missing_genotypes

- **Purpose:** Document a reusable impute missing genotypes capability for Variant analysis workflows. Candidate CAP-0802 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the impute missing genotypes operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.impute_missing_genotypes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `impute_missing_genotypes` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import impute_missing_genotypes

result = impute_missing_genotypes(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0802",
  "operation": "impute_missing_genotypes",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'impute_missing_genotypes', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 45) left_align

- **Purpose:** Document a reusable left align capability for Variant analysis workflows. Candidate CAP-0803 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the left align operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.left_align(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `left_align` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import left_align

result = left_align(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0803",
  "operation": "left_align",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'left_align', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a218...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 46) left_align_mt

- **Purpose:** Document a reusable left align mt capability for Variant analysis workflows. Candidate CAP-0804 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the left align mt operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.left_align_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `left_align_mt` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import left_align_mt

result = left_align_mt(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0804",
  "operation": "left_align_mt",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'left_align_mt', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 47) left_align_variants

- **Purpose:** Document a reusable left align variants capability for Variant analysis workflows. Candidate CAP-0805 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the left align variants operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.left_align_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `left_align_variants` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import left_align_variants

result = left_align_variants(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0805",
  "operation": "left_align_variants",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'left_align_variants', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e0...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 48) merge

- **Purpose:** Document a reusable merge capability for Variant analysis workflows. Candidate CAP-0806 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.merge(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `merge` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import merge

result = merge(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0806",
  "operation": "merge",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'merge', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a21887'},...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 49) merge_mt

- **Purpose:** Document a reusable merge mt capability for Variant analysis workflows. Candidate CAP-0807 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge mt operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.merge_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `merge_mt` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import merge_mt

result = merge_mt(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0807",
  "operation": "merge_mt",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'merge_mt', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a21887...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 50) merge_vcf_files

- **Purpose:** Document a reusable merge vcf files capability for Variant analysis workflows. Candidate CAP-0808 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge vcf files operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.merge_vcf_files(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `merge_vcf_files` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import merge_vcf_files

result = merge_vcf_files(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0808",
  "operation": "merge_vcf_files",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'merge_vcf_files', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a38...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 51) normalize

- **Purpose:** Document a reusable normalize capability for Variant analysis workflows. Candidate CAP-0809 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the normalize operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `generic_operation_name`. Existing-name evidence: `dnakit.standardize.normalize|normalize|vcf_normalize`.

- **API:** `dnakit.variants.normalize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `normalize` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import normalize

result = normalize(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0809",
  "operation": "normalize",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'normalize', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a2188...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 52) normalize_mt

- **Purpose:** Document a reusable normalize mt capability for Variant analysis workflows. Candidate CAP-0810 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the normalize mt operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. Existing-name evidence: `normalize`.

- **API:** `dnakit.variants.normalize_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `normalize_mt` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import normalize_mt

result = normalize_mt(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0810",
  "operation": "normalize_mt",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'normalize_mt', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 53) normalize_vcf

- **Purpose:** Document a reusable normalize vcf capability for Variant analysis workflows. Candidate CAP-0811 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the normalize vcf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. Existing-name evidence: `normalize|vcf_normalize`.

- **API:** `dnakit.variants.normalize_vcf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `normalize_vcf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import normalize_vcf

result = normalize_vcf(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0811",
  "operation": "normalize_vcf",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'normalize_vcf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 54) parse_vcf_records

- **Purpose:** Document a reusable parse vcf records capability for Variant analysis workflows. Candidate CAP-0812 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the parse vcf records operation, and return structured results with provenance and diagnostics. The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.parse_vcf_records(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be VCF text, UTF-8 bytes, a sequence of VCF lines, or a mapping containing `text`; mandatory VCF columns are validated.

- **Output:** The local `parse_vcf_records` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import parse_vcf_records

result = parse_vcf_records(
    input='##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\nchr1\t10\trs1\tA\tG,T\t50\tPASS\tDP=12\tGT\t0/1',
    options={"strict": True}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0812",
  "operation": "parse_vcf_records",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'parse_vcf_records', 'result': {'headers': ['##fileformat=VCFv4.3', '#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\tFORMAT\\ts1'], 'samples': ['s1'], 'records': [{'chrom': 'chr1', 'pos': 10, 'id': 'rs1', 'ref': 'A', 'alt': ['G', 'T'], 'qual': 50.0, 'filter': 'P...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 55) phase_genotypes

- **Purpose:** Document a reusable phase genotypes capability for Variant analysis workflows. Candidate CAP-0813 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the phase genotypes operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.phase_genotypes(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `phase_genotypes` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import phase_genotypes

result = phase_genotypes(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0813",
  "operation": "phase_genotypes",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'phase_genotypes', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a38...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 56) plot_density

- **Purpose:** Document a reusable plot density capability for Variant analysis workflows. Candidate CAP-0814 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the plot density operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.plot_density(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `plot_density` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import plot_density

result = plot_density(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0814",
  "operation": "plot_density",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'plot_density', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 57) plot_density_mt

- **Purpose:** Document a reusable plot density mt capability for Variant analysis workflows. Candidate CAP-0815 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the plot density mt operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.plot_density_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `plot_density_mt` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import plot_density_mt

result = plot_density_mt(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0815",
  "operation": "plot_density_mt",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'plot_density_mt', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a38...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 58) prioritize_candidate_variants

- **Purpose:** Document a reusable prioritize candidate variants capability for Variant analysis workflows. Candidate CAP-0816 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the prioritize candidate variants operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.prioritize_candidate_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `prioritize_candidate_variants` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import prioritize_candidate_variants

result = prioritize_candidate_variants(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0816",
  "operation": "prioritize_candidate_variants",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'prioritize_candidate_variants', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 59) split_multiallelic

- **Purpose:** Document a reusable split multiallelic capability for Variant analysis workflows. Candidate CAP-0817 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split multiallelic operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The dependency-free implementation is available at `src/dnakit/variants/statistics.py` and is exported through the formal DNAKit namespace. It computes descriptive statistics from explicit VCF-like records; it does not call a variant caller, annotation database, reference genome, or imputation backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.split_multiallelic(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input` must be a sequence or mapping containing VCF-like records; each record requires `chrom`, positive 1-based `pos`, `ref`, and `alt`.

- **Output:** The local `split_multiallelic` function returns a completed JSON-serializable envelope with `status`, `operation`, `result`, and explicit provenance. The result is descriptive only: genotypes, annotations, reference alleles, and databases are not inferred.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import split_multiallelic

result = split_multiallelic(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": ["G", "T"], "genotypes": ["0/1", "1/2"]}]},
    options={"remap_genotypes": True}
)
print(result)
```

- **Example output:**

```json
{
  "status": "completed",
  "capability_id": "CAP-0817",
  "operation": "split_multiallelic",
  "result_type": "dict",
  "result_preview": "{'status': 'completed', 'operation': 'split_multiallelic', 'result': {'records': [{'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'G', 'genotypes': ['0/1', '1/0'], 'original_alt_index': 1}, {'chrom': 'chr1', 'pos': 10, 'ref': 'A', 'alt': 'T', 'genotypes': ['0/0', '0/1'], 'original_alt_index': 2}]...",
  "provenance": {
    "implementation": "dnakit.native_variant_statistics",
    "backend": "local_dependency_free",
    "scientific_boundary": "descriptive record-level operation; no caller, annotation database, or imputation backend"
  }
}
```

## 60) split_multiallelic_mt

- **Purpose:** Document a reusable split multiallelic mt capability for Variant analysis workflows. Candidate CAP-0818 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split multiallelic mt operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.split_multiallelic_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `split_multiallelic_mt` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import split_multiallelic_mt

result = split_multiallelic_mt(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0818",
  "operation": "split_multiallelic_mt",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'split_multiallelic_mt', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 61) split_multiallelic_variants

- **Purpose:** Document a reusable split multiallelic variants capability for Variant analysis workflows. Candidate CAP-0819 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split multiallelic variants operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.split_multiallelic_variants(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `split_multiallelic_variants` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import split_multiallelic_variants

result = split_multiallelic_variants(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0819",
  "operation": "split_multiallelic_variants",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'split_multiallelic_variants', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d9...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 62) summarize

- **Purpose:** Document a reusable summarize capability for Variant analysis workflows. Candidate CAP-0820 comes from 14 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the summarize operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.summarize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `summarize` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import summarize

result = summarize(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0820",
  "operation": "summarize",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'summarize', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a2188...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 63) summarize_mt

- **Purpose:** Document a reusable summarize mt capability for Variant analysis workflows. Candidate CAP-0821 comes from 1 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the summarize mt operation, and return structured results with provenance and diagnostics. The optional `variant_type` parameter selects one supported variant (the supported values recorded in the source catalog). The implementation is available at `src/dnakit/variants/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.variants.summarize_mt(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `summarize_mt` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.variants import summarize_mt

result = summarize_mt(
    input={"variants": [{"chrom": "chr1", "pos": 10, "ref": "A", "alt": "G", "genotypes": ["0/1", "1/1", "0/0"]}]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0821",
  "operation": "summarize_mt",
  "result_type": "dict",
  "result_preview": "{'module': 'variants', 'operation': 'summarize_mt', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 1, 'text_length': 104, 'sequence_bases': 0, 'lengths': [], 'sha256': '574c1f57be166487d753e1a2dfee8bc181fbf429595d91ee679e02a383a2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```
