# Input and output formats

This category contains **131** documented candidate capabilities; 131 have a locally callable formal entry (0 partial due to binary-format dependencies; 0 bounded external backends; 112 dependency-adapter only). The documented functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) build_file_manifest

- **Purpose:** Document a reusable build file manifest capability for Input and output formats workflows. Candidate CAP-0308 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the build file manifest operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.build_file_manifest(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `build_file_manifest` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import build_file_manifest

result = build_file_manifest(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0308",
  "operation": "build_file_manifest",
  "result_type": "dict",
  "result_preview": "{'path': None, 'size_bytes': 8, 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232b', 'media_type': None}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 2) chunk_large_file

- **Purpose:** Document a reusable chunk large file capability for Input and output formats workflows. Candidate CAP-0309 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the chunk large file operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.chunk_large_file(input[required], chunk_size[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `chunk_size`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `chunk_large_file` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import chunk_large_file

result = chunk_large_file(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0309",
  "operation": "chunk_large_file",
  "result_type": "list",
  "result_preview": "[b'>a\\nACGT\\n']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 3) clean_temporary_files

- **Purpose:** Document a reusable clean temporary files capability for Input and output formats workflows. Candidate CAP-0310 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the clean temporary files operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.clean_temporary_files(input[required], suffixes[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `suffixes`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `clean_temporary_files` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import clean_temporary_files

result = clean_temporary_files(
    input=Path("nonexistent-temp-directory")
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0310",
  "operation": "clean_temporary_files",
  "result_type": "list",
  "result_preview": "[]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 4) compress

- **Purpose:** Document a reusable compress capability for Input and output formats workflows. Candidate CAP-0311 comes from 25 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compress operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (arrow, bam, bcf, bed, bedgraph, bigbed, bigwig, bim, cram, csv, embl, fam, and 13 more). The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. Existing-name evidence: `compare`.

- **API:** `dnakit.io_formats.compress(input[required], level[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `level`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `compress` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import compress

result = compress(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0311",
  "operation": "compress",
  "result_type": "bytes",
  "result_preview": "<bytes length=28>",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 5) compress_map

- **Purpose:** Document a reusable compress map capability for Input and output formats workflows. Candidate CAP-0312 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compress map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.compress_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compress_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import compress_map

result = compress_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0312",
  "operation": "compress_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'compress_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 6) compress_newick

- **Purpose:** Document a reusable compress newick capability for Input and output formats workflows. Candidate CAP-0313 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compress newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.compress_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compress_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import compress_newick

result = compress_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0313",
  "operation": "compress_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'compress_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 7) compress_nexus

- **Purpose:** Document a reusable compress nexus capability for Input and output formats workflows. Candidate CAP-0314 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compress nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.compress_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compress_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import compress_nexus

result = compress_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0314",
  "operation": "compress_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'compress_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa2239...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 8) compress_paf

- **Purpose:** Document a reusable compress paf capability for Input and output formats workflows. Candidate CAP-0315 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compress paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.compress_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compress_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import compress_paf

result = compress_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0315",
  "operation": "compress_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'compress_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 9) compress_ped

- **Purpose:** Document a reusable compress ped capability for Input and output formats workflows. Candidate CAP-0316 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compress ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.compress_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compress_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import compress_ped

result = compress_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0316",
  "operation": "compress_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'compress_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 10) compress_phylip

- **Purpose:** Document a reusable compress phylip capability for Input and output formats workflows. Candidate CAP-0317 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compress phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.compress_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compress_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import compress_phylip

result = compress_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0317",
  "operation": "compress_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'compress_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 11) compress_twobit

- **Purpose:** Document a reusable compress twobit capability for Input and output formats workflows. Candidate CAP-0318 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compress twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.compress_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compress_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import compress_twobit

result = compress_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0318",
  "operation": "compress_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'compress_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 12) compress_zarr

- **Purpose:** Document a reusable compress zarr capability for Input and output formats workflows. Candidate CAP-0319 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the compress zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.compress_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `compress_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import compress_zarr

result = compress_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0319",
  "operation": "compress_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'compress_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22398...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 13) convert

- **Purpose:** Document a reusable convert capability for Input and output formats workflows. Candidate CAP-0320 comes from 25 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (arrow, bam, bcf, bed, bedgraph, bigbed, bigwig, bim, cram, csv, embl, fam, and 13 more). The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. Existing-name evidence: `id_convert`.

- **API:** `dnakit.io_formats.convert(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `convert` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import convert

result = convert(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0320",
  "operation": "convert",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'convert', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 14) convert_map

- **Purpose:** Document a reusable convert map capability for Input and output formats workflows. Candidate CAP-0321 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.convert_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `convert_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import convert_map

result = convert_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0321",
  "operation": "convert_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'convert_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 15) convert_newick

- **Purpose:** Document a reusable convert newick capability for Input and output formats workflows. Candidate CAP-0322 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.convert_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `convert_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import convert_newick

result = convert_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0322",
  "operation": "convert_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'convert_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa2239...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 16) convert_nexus

- **Purpose:** Document a reusable convert nexus capability for Input and output formats workflows. Candidate CAP-0323 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.convert_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `convert_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import convert_nexus

result = convert_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0323",
  "operation": "convert_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'convert_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22398...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 17) convert_paf

- **Purpose:** Document a reusable convert paf capability for Input and output formats workflows. Candidate CAP-0324 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.convert_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `convert_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import convert_paf

result = convert_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0324",
  "operation": "convert_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'convert_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 18) convert_ped

- **Purpose:** Document a reusable convert ped capability for Input and output formats workflows. Candidate CAP-0325 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.convert_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `convert_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import convert_ped

result = convert_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0325",
  "operation": "convert_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'convert_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 19) convert_phylip

- **Purpose:** Document a reusable convert phylip capability for Input and output formats workflows. Candidate CAP-0326 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.convert_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `convert_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import convert_phylip

result = convert_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0326",
  "operation": "convert_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'convert_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa2239...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 20) convert_twobit

- **Purpose:** Document a reusable convert twobit capability for Input and output formats workflows. Candidate CAP-0327 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.convert_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `convert_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import convert_twobit

result = convert_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0327",
  "operation": "convert_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'convert_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa2239...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 21) convert_zarr

- **Purpose:** Document a reusable convert zarr capability for Input and output formats workflows. Candidate CAP-0328 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the convert zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.convert_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `convert_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import convert_zarr

result = convert_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0328",
  "operation": "convert_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'convert_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 22) decode_metadata_categorical

- **Purpose:** Document a reusable decode metadata categorical capability for Input and output formats workflows. Candidate CAP-0329 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the decode metadata categorical operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.decode_metadata_categorical(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `decode_metadata_categorical` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import decode_metadata_categorical

result = decode_metadata_categorical(
    input={"data": {"group": [0, 1]}, "categories": {"group": ["A", "B"]}}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0329",
  "operation": "decode_metadata_categorical",
  "result_type": "dict",
  "result_preview": "{'group': ['A', 'B']}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 23) decompress

- **Purpose:** Document a reusable decompress capability for Input and output formats workflows. Candidate CAP-0330 comes from 25 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the decompress operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (arrow, bam, bcf, bed, bedgraph, bigbed, bigwig, bim, cram, csv, embl, fam, and 13 more). The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.decompress(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `decompress` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import decompress

result = decompress(
    input=__import__("gzip").compress(b"ACGT")
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0330",
  "operation": "decompress",
  "result_type": "bytes",
  "result_preview": "<bytes length=4>",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 24) decompress_map

- **Purpose:** Document a reusable decompress map capability for Input and output formats workflows. Candidate CAP-0331 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the decompress map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.decompress_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `decompress_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import decompress_map

result = decompress_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0331",
  "operation": "decompress_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'decompress_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa2239...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 25) decompress_newick

- **Purpose:** Document a reusable decompress newick capability for Input and output formats workflows. Candidate CAP-0332 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the decompress newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.decompress_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `decompress_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import decompress_newick

result = decompress_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0332",
  "operation": "decompress_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'decompress_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 26) decompress_nexus

- **Purpose:** Document a reusable decompress nexus capability for Input and output formats workflows. Candidate CAP-0333 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the decompress nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.decompress_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `decompress_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import decompress_nexus

result = decompress_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0333",
  "operation": "decompress_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'decompress_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 27) decompress_paf

- **Purpose:** Document a reusable decompress paf capability for Input and output formats workflows. Candidate CAP-0334 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the decompress paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.decompress_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `decompress_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import decompress_paf

result = decompress_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0334",
  "operation": "decompress_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'decompress_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa2239...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 28) decompress_ped

- **Purpose:** Document a reusable decompress ped capability for Input and output formats workflows. Candidate CAP-0335 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the decompress ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.decompress_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `decompress_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import decompress_ped

result = decompress_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0335",
  "operation": "decompress_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'decompress_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa2239...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 29) decompress_phylip

- **Purpose:** Document a reusable decompress phylip capability for Input and output formats workflows. Candidate CAP-0336 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the decompress phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.decompress_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `decompress_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import decompress_phylip

result = decompress_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0336",
  "operation": "decompress_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'decompress_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 30) decompress_twobit

- **Purpose:** Document a reusable decompress twobit capability for Input and output formats workflows. Candidate CAP-0337 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the decompress twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.decompress_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `decompress_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import decompress_twobit

result = decompress_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0337",
  "operation": "decompress_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'decompress_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 31) decompress_zarr

- **Purpose:** Document a reusable decompress zarr capability for Input and output formats workflows. Candidate CAP-0338 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the decompress zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.decompress_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `decompress_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import decompress_zarr

result = decompress_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0338",
  "operation": "decompress_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'decompress_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 32) detect_compression

- **Purpose:** Document a reusable detect compression capability for Input and output formats workflows. Candidate CAP-0339 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect compression operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.detect_compression(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `detect_compression` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import detect_compression

result = detect_compression(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0339",
  "operation": "detect_compression",
  "result_type": "str",
  "result_preview": "'none'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 33) detect_file_format

- **Purpose:** Document a reusable detect file format capability for Input and output formats workflows. Candidate CAP-0340 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the detect file format operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.detect_file_format(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `detect_file_format` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import detect_file_format

result = detect_file_format(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0340",
  "operation": "detect_file_format",
  "result_type": "str",
  "result_preview": "'fasta'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 34) encode_metadata_categorical

- **Purpose:** Document a reusable encode metadata categorical capability for Input and output formats workflows. Candidate CAP-0342 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the encode metadata categorical operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.encode_metadata_categorical(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `encode_metadata_categorical` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import encode_metadata_categorical

result = encode_metadata_categorical(
    input={"group": ["A", "B", "A"]}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0342",
  "operation": "encode_metadata_categorical",
  "result_type": "dict",
  "result_preview": "{'data': {'group': [0, 1, 0]}, 'categories': {'group': ['A', 'B']}}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 35) estimate_storage_size

- **Purpose:** Document a reusable estimate storage size capability for Input and output formats workflows. Candidate CAP-0343 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the estimate storage size operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.estimate_storage_size(input[required], compression[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `compression`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `estimate_storage_size` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import estimate_storage_size

result = estimate_storage_size(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0343",
  "operation": "estimate_storage_size",
  "result_type": "int",
  "result_preview": "8",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 36) filter

- **Purpose:** Document a reusable filter capability for Input and output formats workflows. Candidate CAP-0344 comes from 25 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (arrow, bam, bcf, bed, bedgraph, bigbed, bigwig, bim, cram, csv, embl, fam, and 13 more). The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.filter(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `filter` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import filter

result = filter(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0344",
  "operation": "filter",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'filter', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232b'...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 37) filter_map

- **Purpose:** Document a reusable filter map capability for Input and output formats workflows. Candidate CAP-0345 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.filter_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `filter_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import filter_map

result = filter_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0345",
  "operation": "filter_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'filter_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 38) filter_newick

- **Purpose:** Document a reusable filter newick capability for Input and output formats workflows. Candidate CAP-0346 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.filter_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `filter_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import filter_newick

result = filter_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0346",
  "operation": "filter_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'filter_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22398...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 39) filter_nexus

- **Purpose:** Document a reusable filter nexus capability for Input and output formats workflows. Candidate CAP-0347 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.filter_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `filter_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import filter_nexus

result = filter_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0347",
  "operation": "filter_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'filter_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 40) filter_paf

- **Purpose:** Document a reusable filter paf capability for Input and output formats workflows. Candidate CAP-0348 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.filter_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `filter_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import filter_paf

result = filter_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0348",
  "operation": "filter_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'filter_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 41) filter_ped

- **Purpose:** Document a reusable filter ped capability for Input and output formats workflows. Candidate CAP-0349 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.filter_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `filter_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import filter_ped

result = filter_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0349",
  "operation": "filter_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'filter_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 42) filter_phylip

- **Purpose:** Document a reusable filter phylip capability for Input and output formats workflows. Candidate CAP-0350 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.filter_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `filter_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import filter_phylip

result = filter_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0350",
  "operation": "filter_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'filter_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22398...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 43) filter_twobit

- **Purpose:** Document a reusable filter twobit capability for Input and output formats workflows. Candidate CAP-0351 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.filter_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `filter_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import filter_twobit

result = filter_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0351",
  "operation": "filter_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'filter_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22398...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 44) filter_zarr

- **Purpose:** Document a reusable filter zarr capability for Input and output formats workflows. Candidate CAP-0352 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the filter zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.filter_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `filter_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import filter_zarr

result = filter_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0352",
  "operation": "filter_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'filter_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 45) get_file_checksum

- **Purpose:** Document a reusable get file checksum capability for Input and output formats workflows. Candidate CAP-0353 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the get file checksum operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.get_file_checksum(input[required], algorithm[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `algorithm`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `get_file_checksum` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import get_file_checksum

result = get_file_checksum(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0353",
  "operation": "get_file_checksum",
  "result_type": "str",
  "result_preview": "'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232b'",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 46) index

- **Purpose:** Document a reusable index capability for Input and output formats workflows. Candidate CAP-0354 comes from 25 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the index operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (arrow, bam, bcf, bed, bedgraph, bigbed, bigwig, bim, cram, csv, embl, fam, and 13 more). The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.index(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `index` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import index

result = index(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0354",
  "operation": "index",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'index', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232b'}...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 47) index_map

- **Purpose:** Document a reusable index map capability for Input and output formats workflows. Candidate CAP-0355 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the index map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.index_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `index_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import index_map

result = index_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0355",
  "operation": "index_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'index_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 48) index_newick

- **Purpose:** Document a reusable index newick capability for Input and output formats workflows. Candidate CAP-0356 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the index newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.index_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `index_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import index_newick

result = index_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0356",
  "operation": "index_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'index_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 49) index_nexus

- **Purpose:** Document a reusable index nexus capability for Input and output formats workflows. Candidate CAP-0357 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the index nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.index_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `index_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import index_nexus

result = index_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0357",
  "operation": "index_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'index_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 50) index_paf

- **Purpose:** Document a reusable index paf capability for Input and output formats workflows. Candidate CAP-0358 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the index paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.index_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `index_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import index_paf

result = index_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0358",
  "operation": "index_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'index_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 51) index_ped

- **Purpose:** Document a reusable index ped capability for Input and output formats workflows. Candidate CAP-0359 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the index ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.index_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `index_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import index_ped

result = index_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0359",
  "operation": "index_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'index_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 52) index_phylip

- **Purpose:** Document a reusable index phylip capability for Input and output formats workflows. Candidate CAP-0360 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the index phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.index_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `index_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import index_phylip

result = index_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0360",
  "operation": "index_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'index_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 53) index_twobit

- **Purpose:** Document a reusable index twobit capability for Input and output formats workflows. Candidate CAP-0361 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the index twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.index_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `index_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import index_twobit

result = index_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0361",
  "operation": "index_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'index_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 54) index_zarr

- **Purpose:** Document a reusable index zarr capability for Input and output formats workflows. Candidate CAP-0362 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the index zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.index_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `index_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import index_zarr

result = index_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0362",
  "operation": "index_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'index_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 55) list_file_records

- **Purpose:** Document a reusable list file records capability for Input and output formats workflows. Candidate CAP-0363 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the list file records operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.list_file_records(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `list_file_records` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import list_file_records

result = list_file_records(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0363",
  "operation": "list_file_records",
  "result_type": "list",
  "result_preview": "['>a', 'ACGT']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 56) merge

- **Purpose:** Document a reusable merge capability for Input and output formats workflows. Candidate CAP-0364 comes from 25 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (arrow, bam, bcf, bed, bedgraph, bigbed, bigwig, bim, cram, csv, embl, fam, and 13 more). The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.merge(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `merge` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import merge

result = merge(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0364",
  "operation": "merge",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'merge', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232b'}...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 57) merge_map

- **Purpose:** Document a reusable merge map capability for Input and output formats workflows. Candidate CAP-0365 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.merge_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `merge_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import merge_map

result = merge_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0365",
  "operation": "merge_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'merge_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 58) merge_metadata_tables

- **Purpose:** Document a reusable merge metadata tables capability for Input and output formats workflows. Candidate CAP-0366 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge metadata tables operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. Existing-name evidence: `merge_metadata`.

- **API:** `dnakit.io_formats.merge_metadata_tables(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `merge_metadata_tables` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import merge_metadata_tables

result = merge_metadata_tables(
    input=[[{"sample_id": "sample1"}], [{"sample_id": "sample2"}]]
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0366",
  "operation": "merge_metadata_tables",
  "result_type": "list",
  "result_preview": "[{'sample_id': 'sample1'}, {'sample_id': 'sample2'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 59) merge_newick

- **Purpose:** Document a reusable merge newick capability for Input and output formats workflows. Candidate CAP-0367 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.merge_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `merge_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import merge_newick

result = merge_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0367",
  "operation": "merge_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'merge_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 60) merge_nexus

- **Purpose:** Document a reusable merge nexus capability for Input and output formats workflows. Candidate CAP-0368 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.merge_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `merge_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import merge_nexus

result = merge_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0368",
  "operation": "merge_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'merge_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 61) merge_paf

- **Purpose:** Document a reusable merge paf capability for Input and output formats workflows. Candidate CAP-0369 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.merge_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `merge_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import merge_paf

result = merge_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0369",
  "operation": "merge_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'merge_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 62) merge_ped

- **Purpose:** Document a reusable merge ped capability for Input and output formats workflows. Candidate CAP-0370 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.merge_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `merge_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import merge_ped

result = merge_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0370",
  "operation": "merge_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'merge_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 63) merge_phylip

- **Purpose:** Document a reusable merge phylip capability for Input and output formats workflows. Candidate CAP-0371 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.merge_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `merge_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import merge_phylip

result = merge_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0371",
  "operation": "merge_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'merge_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 64) merge_twobit

- **Purpose:** Document a reusable merge twobit capability for Input and output formats workflows. Candidate CAP-0372 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.merge_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `merge_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import merge_twobit

result = merge_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0372",
  "operation": "merge_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'merge_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 65) merge_zarr

- **Purpose:** Document a reusable merge zarr capability for Input and output formats workflows. Candidate CAP-0373 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the merge zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.merge_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `merge_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import merge_zarr

result = merge_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0373",
  "operation": "merge_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'merge_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 66) normalize_sample_metadata

- **Purpose:** Document a reusable normalize sample metadata capability for Input and output formats workflows. Candidate CAP-0374 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the normalize sample metadata operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.normalize_sample_metadata(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `normalize_sample_metadata` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import normalize_sample_metadata

result = normalize_sample_metadata(
    input=[{"sample_id": "sample1", "group": "A"}]
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0374",
  "operation": "normalize_sample_metadata",
  "result_type": "list",
  "result_preview": "[{'sample_id': 'sample1', 'group': 'A'}]",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 67) read_map

- **Purpose:** Document a reusable read map capability for Input and output formats workflows. Candidate CAP-0376 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the read map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. Existing-name evidence: `read_agp`.

- **API:** `dnakit.io_formats.read_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `read_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import read_map

result = read_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0376",
  "operation": "read_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'read_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 68) read_newick

- **Purpose:** Document a reusable read newick capability for Input and output formats workflows. Candidate CAP-0377 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the read newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.read_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `read_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import read_newick

result = read_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0377",
  "operation": "read_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'read_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 69) read_nexus

- **Purpose:** Document a reusable read nexus capability for Input and output formats workflows. Candidate CAP-0378 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the read nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.read_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `read_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import read_nexus

result = read_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0378",
  "operation": "read_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'read_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 70) read_paf

- **Purpose:** Document a reusable read paf capability for Input and output formats workflows. Candidate CAP-0379 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the read paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.read_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `read_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import read_paf

result = read_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0379",
  "operation": "read_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'read_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 71) read_ped

- **Purpose:** Document a reusable read ped capability for Input and output formats workflows. Candidate CAP-0380 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the read ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. Existing-name evidence: `read_bed`.

- **API:** `dnakit.io_formats.read_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `read_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import read_ped

result = read_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0380",
  "operation": "read_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'read_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 72) read_phylip

- **Purpose:** Document a reusable read phylip capability for Input and output formats workflows. Candidate CAP-0381 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the read phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.read_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `read_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import read_phylip

result = read_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0381",
  "operation": "read_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'read_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 73) read_twobit

- **Purpose:** Document a reusable read twobit capability for Input and output formats workflows. Candidate CAP-0382 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the read twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.read_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `read_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import read_twobit

result = read_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0382",
  "operation": "read_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'read_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 74) read_zarr

- **Purpose:** Document a reusable read zarr capability for Input and output formats workflows. Candidate CAP-0383 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the read zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.read_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `read_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import read_zarr

result = read_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0383",
  "operation": "read_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'read_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 75) repack_compressed_file

- **Purpose:** Document a reusable repack compressed file capability for Input and output formats workflows. Candidate CAP-0384 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the repack compressed file operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.repack_compressed_file(input[required], level[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `level`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `repack_compressed_file` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import repack_compressed_file

result = repack_compressed_file(
    input=__import__("gzip").compress(b"ACGT")
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0384",
  "operation": "repack_compressed_file",
  "result_type": "bytes",
  "result_preview": "<bytes length=24>",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 76) resume_download

- **Purpose:** Document a reusable resume download capability for Input and output formats workflows. Candidate CAP-0385 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the resume download operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.resume_download(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `resume_download` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import resume_download

result = resume_download(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0385",
  "operation": "resume_download",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'resume_download', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 77) sort

- **Purpose:** Document a reusable sort capability for Input and output formats workflows. Candidate CAP-0386 comes from 25 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sort operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (arrow, bam, bcf, bed, bedgraph, bigbed, bigwig, bim, cram, csv, embl, fam, and 13 more). The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.sort(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `sort` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import sort

result = sort(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0386",
  "operation": "sort",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'sort', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232b'},...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 78) sort_map

- **Purpose:** Document a reusable sort map capability for Input and output formats workflows. Candidate CAP-0387 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sort map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.sort_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `sort_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import sort_map

result = sort_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0387",
  "operation": "sort_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'sort_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 79) sort_newick

- **Purpose:** Document a reusable sort newick capability for Input and output formats workflows. Candidate CAP-0388 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sort newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.sort_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `sort_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import sort_newick

result = sort_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0388",
  "operation": "sort_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'sort_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 80) sort_nexus

- **Purpose:** Document a reusable sort nexus capability for Input and output formats workflows. Candidate CAP-0389 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sort nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.sort_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `sort_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import sort_nexus

result = sort_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0389",
  "operation": "sort_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'sort_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 81) sort_paf

- **Purpose:** Document a reusable sort paf capability for Input and output formats workflows. Candidate CAP-0390 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sort paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.sort_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `sort_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import sort_paf

result = sort_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0390",
  "operation": "sort_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'sort_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 82) sort_ped

- **Purpose:** Document a reusable sort ped capability for Input and output formats workflows. Candidate CAP-0391 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sort ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.sort_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `sort_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import sort_ped

result = sort_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0391",
  "operation": "sort_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'sort_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 83) sort_phylip

- **Purpose:** Document a reusable sort phylip capability for Input and output formats workflows. Candidate CAP-0392 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sort phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.sort_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `sort_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import sort_phylip

result = sort_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0392",
  "operation": "sort_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'sort_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 84) sort_twobit

- **Purpose:** Document a reusable sort twobit capability for Input and output formats workflows. Candidate CAP-0393 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sort twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.sort_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `sort_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import sort_twobit

result = sort_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0393",
  "operation": "sort_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'sort_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 85) sort_zarr

- **Purpose:** Document a reusable sort zarr capability for Input and output formats workflows. Candidate CAP-0394 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the sort zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.sort_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `sort_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import sort_zarr

result = sort_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0394",
  "operation": "sort_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'sort_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 86) split_map

- **Purpose:** Document a reusable split map capability for Input and output formats workflows. Candidate CAP-0396 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.split_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `split_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import split_map

result = split_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0396",
  "operation": "split_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'split_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 87) split_newick

- **Purpose:** Document a reusable split newick capability for Input and output formats workflows. Candidate CAP-0397 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.split_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `split_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import split_newick

result = split_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0397",
  "operation": "split_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'split_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 88) split_nexus

- **Purpose:** Document a reusable split nexus capability for Input and output formats workflows. Candidate CAP-0398 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.split_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `split_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import split_nexus

result = split_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0398",
  "operation": "split_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'split_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 89) split_paf

- **Purpose:** Document a reusable split paf capability for Input and output formats workflows. Candidate CAP-0399 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.split_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `split_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import split_paf

result = split_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0399",
  "operation": "split_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'split_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 90) split_ped

- **Purpose:** Document a reusable split ped capability for Input and output formats workflows. Candidate CAP-0400 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.split_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `split_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import split_ped

result = split_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0400",
  "operation": "split_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'split_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 91) split_phylip

- **Purpose:** Document a reusable split phylip capability for Input and output formats workflows. Candidate CAP-0401 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.split_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `split_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import split_phylip

result = split_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0401",
  "operation": "split_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'split_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 92) split_twobit

- **Purpose:** Document a reusable split twobit capability for Input and output formats workflows. Candidate CAP-0402 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.split_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `split_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import split_twobit

result = split_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0402",
  "operation": "split_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'split_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 93) split_zarr

- **Purpose:** Document a reusable split zarr capability for Input and output formats workflows. Candidate CAP-0403 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the split zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.split_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `split_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import split_zarr

result = split_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0403",
  "operation": "split_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'split_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 94) stream

- **Purpose:** Document a reusable stream capability for Input and output formats workflows. Candidate CAP-0404 comes from 25 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the stream operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (arrow, bam, bcf, bed, bedgraph, bigbed, bigwig, bim, cram, csv, embl, fam, and 13 more). The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.stream(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `stream` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import stream

result = stream(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0404",
  "operation": "stream",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'stream', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b232b'...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 95) stream_map

- **Purpose:** Document a reusable stream map capability for Input and output formats workflows. Candidate CAP-0405 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the stream map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.stream_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `stream_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import stream_map

result = stream_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0405",
  "operation": "stream_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'stream_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 96) stream_newick

- **Purpose:** Document a reusable stream newick capability for Input and output formats workflows. Candidate CAP-0406 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the stream newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.stream_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `stream_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import stream_newick

result = stream_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0406",
  "operation": "stream_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'stream_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22398...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 97) stream_nexus

- **Purpose:** Document a reusable stream nexus capability for Input and output formats workflows. Candidate CAP-0407 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the stream nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.stream_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `stream_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import stream_nexus

result = stream_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0407",
  "operation": "stream_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'stream_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 98) stream_paf

- **Purpose:** Document a reusable stream paf capability for Input and output formats workflows. Candidate CAP-0408 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the stream paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.stream_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `stream_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import stream_paf

result = stream_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0408",
  "operation": "stream_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'stream_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 99) stream_ped

- **Purpose:** Document a reusable stream ped capability for Input and output formats workflows. Candidate CAP-0409 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the stream ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.stream_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `stream_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import stream_ped

result = stream_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0409",
  "operation": "stream_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'stream_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 100) stream_phylip

- **Purpose:** Document a reusable stream phylip capability for Input and output formats workflows. Candidate CAP-0410 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the stream phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.stream_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `stream_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import stream_phylip

result = stream_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0410",
  "operation": "stream_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'stream_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22398...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 101) stream_records

- **Purpose:** Document a reusable stream records capability for Input and output formats workflows. Candidate CAP-0411 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the stream records operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. Existing-name evidence: `trim_record`.

- **API:** `dnakit.io_formats.stream_records(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `stream_records` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import stream_records

result = stream_records(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0411",
  "operation": "stream_records",
  "result_type": "list",
  "result_preview": "['>a', 'ACGT']",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 102) stream_twobit

- **Purpose:** Document a reusable stream twobit capability for Input and output formats workflows. Candidate CAP-0412 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the stream twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.stream_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `stream_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import stream_twobit

result = stream_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0412",
  "operation": "stream_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'stream_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22398...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 103) stream_zarr

- **Purpose:** Document a reusable stream zarr capability for Input and output formats workflows. Candidate CAP-0413 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the stream zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.stream_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `stream_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import stream_zarr

result = stream_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0413",
  "operation": "stream_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'stream_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 104) summarize

- **Purpose:** Document a reusable summarize capability for Input and output formats workflows. Candidate CAP-0414 comes from 25 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the summarize operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (arrow, bam, bcf, bed, bedgraph, bigbed, bigwig, bim, cram, csv, embl, fam, and 13 more). The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `merged_dimension`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.summarize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `summarize` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import summarize

result = summarize(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0414",
  "operation": "summarize",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'summarize', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 105) summarize_map

- **Purpose:** Document a reusable summarize map capability for Input and output formats workflows. Candidate CAP-0415 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the summarize map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.summarize_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `summarize_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import summarize_map

result = summarize_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0415",
  "operation": "summarize_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'summarize_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22398...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 106) summarize_newick

- **Purpose:** Document a reusable summarize newick capability for Input and output formats workflows. Candidate CAP-0416 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the summarize newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.summarize_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `summarize_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import summarize_newick

result = summarize_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0416",
  "operation": "summarize_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'summarize_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 107) summarize_nexus

- **Purpose:** Document a reusable summarize nexus capability for Input and output formats workflows. Candidate CAP-0417 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the summarize nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.summarize_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `summarize_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import summarize_nexus

result = summarize_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0417",
  "operation": "summarize_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'summarize_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 108) summarize_paf

- **Purpose:** Document a reusable summarize paf capability for Input and output formats workflows. Candidate CAP-0418 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the summarize paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.summarize_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `summarize_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import summarize_paf

result = summarize_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0418",
  "operation": "summarize_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'summarize_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22398...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 109) summarize_ped

- **Purpose:** Document a reusable summarize ped capability for Input and output formats workflows. Candidate CAP-0419 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the summarize ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.summarize_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `summarize_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import summarize_ped

result = summarize_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0419",
  "operation": "summarize_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'summarize_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22398...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 110) summarize_phylip

- **Purpose:** Document a reusable summarize phylip capability for Input and output formats workflows. Candidate CAP-0420 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the summarize phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.summarize_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `summarize_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import summarize_phylip

result = summarize_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0420",
  "operation": "summarize_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'summarize_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 111) summarize_twobit

- **Purpose:** Document a reusable summarize twobit capability for Input and output formats workflows. Candidate CAP-0421 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the summarize twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.summarize_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `summarize_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import summarize_twobit

result = summarize_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0421",
  "operation": "summarize_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'summarize_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 112) summarize_zarr

- **Purpose:** Document a reusable summarize zarr capability for Input and output formats workflows. Candidate CAP-0422 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the summarize zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.summarize_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `summarize_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import summarize_zarr

result = summarize_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0422",
  "operation": "summarize_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'summarize_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa2239...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 113) validate

- **Purpose:** Document a reusable validate capability for Input and output formats workflows. Candidate CAP-0423 comes from 25 source row(s) and is currently marked `review_needed`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate operation, and return structured results with provenance and diagnostics. The optional `format` parameter selects one supported variant (arrow, bam, bcf, bed, bedgraph, bigbed, bigwig, bim, cram, csv, embl, fam, and 13 more). The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `generic_operation_name`. Existing-name evidence: `dnakit.standardize.validate|validate|validate_set`.

- **API:** `dnakit.io_formats.validate(input[required])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: no. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `validate` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import validate

result = validate(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0423",
  "operation": "validate",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'size_bytes': 8, 'line_count': 2, 'error': None}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 114) validate_map

- **Purpose:** Document a reusable validate map capability for Input and output formats workflows. Candidate CAP-0424 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. Existing-name evidence: `validate`.

- **API:** `dnakit.io_formats.validate_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `validate_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import validate_map

result = validate_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0424",
  "operation": "validate_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'validate_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 115) validate_metadata_schema

- **Purpose:** Document a reusable validate metadata schema capability for Input and output formats workflows. Candidate CAP-0425 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate metadata schema operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. Existing-name evidence: `validate_metadata`.

- **API:** `dnakit.io_formats.validate_metadata_schema(input[required], required[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `required`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `validate_metadata_schema` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import validate_metadata_schema

result = validate_metadata_schema(
    input={"sample_id": "sample1"}
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0425",
  "operation": "validate_metadata_schema",
  "result_type": "dict",
  "result_preview": "{'valid': True, 'missing': [], 'keys': ['sample_id']}",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 116) validate_newick

- **Purpose:** Document a reusable validate newick capability for Input and output formats workflows. Candidate CAP-0426 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.validate_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `validate_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import validate_newick

result = validate_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0426",
  "operation": "validate_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'validate_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 117) validate_nexus

- **Purpose:** Document a reusable validate nexus capability for Input and output formats workflows. Candidate CAP-0427 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.validate_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `validate_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import validate_nexus

result = validate_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0427",
  "operation": "validate_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'validate_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa2239...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 118) validate_paf

- **Purpose:** Document a reusable validate paf capability for Input and output formats workflows. Candidate CAP-0428 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. Existing-name evidence: `validate`.

- **API:** `dnakit.io_formats.validate_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `validate_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import validate_paf

result = validate_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0428",
  "operation": "validate_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'validate_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 119) validate_ped

- **Purpose:** Document a reusable validate ped capability for Input and output formats workflows. Candidate CAP-0429 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. Existing-name evidence: `validate|validate_set`.

- **API:** `dnakit.io_formats.validate_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `validate_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import validate_ped

result = validate_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0429",
  "operation": "validate_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'validate_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 120) validate_phylip

- **Purpose:** Document a reusable validate phylip capability for Input and output formats workflows. Candidate CAP-0430 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.validate_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `validate_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import validate_phylip

result = validate_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0430",
  "operation": "validate_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'validate_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 121) validate_twobit

- **Purpose:** Document a reusable validate twobit capability for Input and output formats workflows. Candidate CAP-0431 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.validate_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `validate_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import validate_twobit

result = validate_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0431",
  "operation": "validate_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'validate_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 122) validate_zarr

- **Purpose:** Document a reusable validate zarr capability for Input and output formats workflows. Candidate CAP-0432 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the validate zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.validate_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `validate_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import validate_zarr

result = validate_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0432",
  "operation": "validate_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'validate_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa22398...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 123) verify_checksum

- **Purpose:** Document a reusable verify checksum capability for Input and output formats workflows. Candidate CAP-0433 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the verify checksum operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.verify_checksum(input[required], expected[required], algorithm[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`, `expected`; optional arguments: `algorithm`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local `verify_checksum` function returns its native Python value (for example, a list, mapping, string, bytes, or record object) rather than a formal DNAKit result envelope. The concrete example output below shows the observed type and a bounded preview.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import verify_checksum

result = verify_checksum(
    input='>a\nACGT\n',
    expected=__import__("hashlib").sha256(b">a\nACGT\n").hexdigest()
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_implemented",
  "capability_id": "CAP-0433",
  "operation": "verify_checksum",
  "result_type": "bool",
  "result_preview": "True",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_local"
  }
}
```

## 124) write_map

- **Purpose:** Document a reusable write map capability for Input and output formats workflows. Candidate CAP-0435 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the write map operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. Existing-name evidence: `write_agp`.

- **API:** `dnakit.io_formats.write_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `write_map` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import write_map

result = write_map(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0435",
  "operation": "write_map",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'write_map', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 125) write_newick

- **Purpose:** Document a reusable write newick capability for Input and output formats workflows. Candidate CAP-0436 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the write newick operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.write_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `write_newick` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import write_newick

result = write_newick(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0436",
  "operation": "write_newick",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'write_newick', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 126) write_nexus

- **Purpose:** Document a reusable write nexus capability for Input and output formats workflows. Candidate CAP-0437 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the write nexus operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.write_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `write_nexus` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import write_nexus

result = write_nexus(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0437",
  "operation": "write_nexus",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'write_nexus', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 127) write_paf

- **Purpose:** Document a reusable write paf capability for Input and output formats workflows. Candidate CAP-0438 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the write paf operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.write_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `write_paf` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import write_paf

result = write_paf(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0438",
  "operation": "write_paf",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'write_paf', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 128) write_ped

- **Purpose:** Document a reusable write ped capability for Input and output formats workflows. Candidate CAP-0439 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the write ped operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. Existing-name evidence: `write_bed`.

- **API:** `dnakit.io_formats.write_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `write_ped` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import write_ped

result = write_ped(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0439",
  "operation": "write_ped",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'write_ped', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b23...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 129) write_phylip

- **Purpose:** Document a reusable write phylip capability for Input and output formats workflows. Candidate CAP-0440 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the write phylip operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.write_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `write_phylip` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import write_phylip

result = write_phylip(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0440",
  "operation": "write_phylip",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'write_phylip', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 130) write_twobit

- **Purpose:** Document a reusable write twobit capability for Input and output formats workflows. Candidate CAP-0441 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the write twobit operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.write_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `write_twobit` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import write_twobit

result = write_twobit(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0441",
  "operation": "write_twobit",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'write_twobit', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```

## 131) write_zarr

- **Purpose:** Document a reusable write zarr capability for Input and output formats workflows. Candidate CAP-0442 comes from 1 source row(s) and is currently marked `new_candidate`; the candidate keeps input-type, backend, format, platform, provider, assay, or biological-entity differences inside one typed interface rather than creating one function per variant.

- **Method:** Validate the typed input envelope, normalize records and metadata, perform the write zarr operation, and return structured results with provenance and diagnostics. The implementation is available at `src/dnakit/io_formats/__init__.py` and is exported through the formal DNAKit namespace. Agent/MCP exposure is available when the function's JSON call plan is compatible. This is a local adapter: it validates the typed input and reports the backend required for scientific execution; it does not claim to replace that backend. Review basis: `no_exact_public_match`. No exact public-function name was found in the comparison.

- **API:** `dnakit.io_formats.write_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Required arguments: `input`; optional arguments: `options`. IO arguments accept text, bytes, mappings, paths, or materialized records as required by the signature.

- **Output:** The local adapter `write_zarr` returns a JSON-serializable envelope containing `status`, `module`, `operation`, `input_summary`, `provenance`, and either a local result or `required_dependency`. It does not fabricate a backend scientific result. The adapter returns a structured dependency status when its scientific backend is unavailable.

- **Example code:**

```python
from pathlib import Path

from dnakit.io_formats import write_zarr

result = write_zarr(
    input='>a\nACGT\n'
)
print(result)
```

- **Example output:**

```json
{
  "status": "candidate_adapter",
  "capability_id": "CAP-0442",
  "operation": "write_zarr",
  "result_type": "dict",
  "result_preview": "{'module': 'io_formats', 'operation': 'write_zarr', 'provenance': {'implementation': 'local_candidate_adapter', 'backend': 'local'}, 'input_summary': {'record_count': 2, 'text_length': 8, 'sequence_bases': 4, 'lengths': [4], 'sha256': 'ec93753459551ffadb17f616d5e6bd45de641c045083fac84c4fa223987b2...",
  "provenance": {
    "schema": "candidate_registry",
    "implementation": "formal_adapter"
  }
}
```
