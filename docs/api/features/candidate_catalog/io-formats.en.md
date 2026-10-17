# Input and output formats

This category contains **131** documented capabilities, all callable through the formal DNAKit namespace. The 96 functions described in the format sections compute real results and none of them returns `status: "not_implemented"` unconditionally: the 88 `<operation>_<format>` functions (11 operations `read`, `write`, `validate`, `stream`, `summarize`, `index`, `merge`, `split`, `filter`, `sort`, `convert` for each of 8 formats: PLINK `map` and `ped`, `paf`, `newick`, `nexus`, `phylip`, UCSC `twobit` and `zarr`), the 7 generic functions `convert`, `filter`, `index`, `merge`, `sort`, `stream` and `summarize` that detect the format and call the matching one, and `resume_download`. The 11 `*_zarr` functions need the optional `zarr>=3` and `numpy` packages and return `status: "not_implemented"` (with `required_dependency` and `message`) only while these are not installed; `resume_download` needs network access to the server it downloads from. The format functions write nothing unless `output`/`output_dir` is given and are bounded by `max_bytes`/`max_records`/`max_cells`/`max_elements` limits. The remaining 35 functions (compression, checksums, manifests, metadata tables, ...) are described in their own sections. The functions are exported through the formal DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

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

- **Purpose:** Compress a PLINK `.map` marker table (text, bytes or a file) with deterministic gzip and return it base64-encoded.

- **Method:** The input is turned into bytes (see Input) and compressed with gzip at `options['level']` (1-9, default 6). The gzip header is fixed (modification time 0, no file name), so the same input and level always give the same `data`. The result is base64-encoded and returned with the sizes and the SHA-256 of the uncompressed bytes. The content is never parsed or validated as a PLINK `.map` marker table, and binary content is kept bit-exact; `decompress_map` restores the exact bytes. Nothing is written unless `options['output']` is given: then the raw gzip bytes are written to that file (an existing file needs `overwrite=True`).

- **API:** `dnakit.io_formats.compress_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: text (`str`, encoded as UTF-8), `bytes`/`bytearray`, or the path (`str` or `Path`) of an existing file whose bytes are compressed (a directory raises `InputFormatError`; only the zarr functions accept one). A one-line string that is an existing path is read as a path. A `Path`, or a string with a known file suffix such as `.2bit`, that does not exist raises `InputFormatError` (code `IO_FORMATS_FILE_NOT_FOUND`) instead of being compressed as text. Mappings, lists and numbers are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `level`: int 1-9, default 6: gzip compression level
    - `output`: path of a .gz file to write (optional; nothing is written otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `encoding` (`"gzip+base64"`), `data` (the base64 text), `original_bytes` (size before compression), `compressed_bytes` (size of the gzip stream), `sha256_original` (SHA-256 of the uncompressed bytes) and `output_path` when `options['output']` is given. Unusable inputs raise `InputFormatError`; invalid or unknown options (for example `level` outside 1-9), an existing output file without `overwrite=True` and inputs above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import compress_map

envelope = compress_map('1\trs1\t0\t1000\n1\trs2\t0\t2000\n', verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAAzPkLCo25DTgNDQwMOAyBHKMgBwjEAcAAo9oVBoAAAA=",
    "original_bytes": 26,
    "compressed_bytes": 41,
    "sha256_original": "0d63a7da7a1761ff406f984298d593a88a953100925de4daf0f0b345438301e9"
  }
}
```

## 6) compress_newick

- **Purpose:** Compress a Newick tree file (text, bytes or a file) with deterministic gzip and return it base64-encoded.

- **Method:** The input is turned into bytes (see Input) and compressed with gzip at `options['level']` (1-9, default 6). The gzip header is fixed (modification time 0, no file name), so the same input and level always give the same `data`. The result is base64-encoded and returned with the sizes and the SHA-256 of the uncompressed bytes. The content is never parsed or validated as a Newick tree file, and binary content is kept bit-exact; `decompress_newick` restores the exact bytes. Nothing is written unless `options['output']` is given: then the raw gzip bytes are written to that file (an existing file needs `overwrite=True`).

- **API:** `dnakit.io_formats.compress_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: text (`str`, encoded as UTF-8), `bytes`/`bytearray`, or the path (`str` or `Path`) of an existing file whose bytes are compressed (a directory raises `InputFormatError`; only the zarr functions accept one). A one-line string that is an existing path is read as a path. A `Path`, or a string with a known file suffix such as `.2bit`, that does not exist raises `InputFormatError` (code `IO_FORMATS_FILE_NOT_FOUND`) instead of being compressed as text. Mappings, lists and numbers are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `level`: int 1-9, default 6: gzip compression level
    - `output`: path of a .gz file to write (optional; nothing is written otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `encoding` (`"gzip+base64"`), `data` (the base64 text), `original_bytes` (size before compression), `compressed_bytes` (size of the gzip stream), `sha256_original` (SHA-256 of the uncompressed bytes) and `output_path` when `options['output']` is given. Unusable inputs raise `InputFormatError`; invalid or unknown options (for example `level` outside 1-9), an existing output file without `overwrite=True` and inputs above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import compress_newick

envelope = compress_newick('((A:1,B:2):1,C:3);', verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAA9PQcLQy1HGyMtIEUs5WxprWADsnB7wSAAAA",
    "original_bytes": 18,
    "compressed_bytes": 36,
    "sha256_original": "1687864f64fa679bf85d6ff2ba780918d7c19e3965bcf14f658c9ee9c568c968"
  }
}
```

## 7) compress_nexus

- **Purpose:** Compress a NEXUS file (text, bytes or a file) with deterministic gzip and return it base64-encoded.

- **Method:** The input is turned into bytes (see Input) and compressed with gzip at `options['level']` (1-9, default 6). The gzip header is fixed (modification time 0, no file name), so the same input and level always give the same `data`. The result is base64-encoded and returned with the sizes and the SHA-256 of the uncompressed bytes. The content is never parsed or validated as a NEXUS file, and binary content is kept bit-exact; `decompress_nexus` restores the exact bytes. Nothing is written unless `options['output']` is given: then the raw gzip bytes are written to that file (an existing file needs `overwrite=True`).

- **API:** `dnakit.io_formats.compress_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: text (`str`, encoded as UTF-8), `bytes`/`bytearray`, or the path (`str` or `Path`) of an existing file whose bytes are compressed (a directory raises `InputFormatError`; only the zarr functions accept one). A one-line string that is an existing path is read as a path. A `Path`, or a string with a known file suffix such as `.2bit`, that does not exist raises `InputFormatError` (code `IO_FORMATS_FILE_NOT_FOUND`) instead of being compressed as text. Mappings, lists and numbers are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `level`: int 1-9, default 6: gzip compression level
    - `output`: path of a .gz file to write (optional; nothing is written otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `encoding` (`"gzip+base64"`), `data` (the base64 text), `original_bytes` (size before compression), `compressed_bytes` (size of the gzip stream), `sha256_original` (SHA-256 of the uncompressed bytes) and `output_path` when `options['output']` is given. Unusable inputs raise `InputFormatError`; invalid or unknown options (for example `level` outside 1-9), an existing output file without `overwrite=True` and inputs above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import compress_nexus

envelope = compress_nexus('#NEXUS\nBEGIN TAXA;\n  DIMENSIONS NTAX=2;\n  TAXLABELS A B;\nEND;\n', verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAA1P2c40IDeZycnX39FMIcYxwtOZSUHDx9HX1C/b09wtW8AOK2RqBBIEMH0cnV59gBUcFJ2suVz8Xay4AkYWK6z4AAAA=",
    "original_bytes": 62,
    "compressed_bytes": 77,
    "sha256_original": "80a4f6422311b951751d59815bb5a36465fe705e84fb840ecb8f30d8c4a36a8e"
  }
}
```

## 8) compress_paf

- **Purpose:** Compress a PAF alignment table (text, bytes or a file) with deterministic gzip and return it base64-encoded.

- **Method:** The input is turned into bytes (see Input) and compressed with gzip at `options['level']` (1-9, default 6). The gzip header is fixed (modification time 0, no file name), so the same input and level always give the same `data`. The result is base64-encoded and returned with the sizes and the SHA-256 of the uncompressed bytes. The content is never parsed or validated as a PAF alignment table, and binary content is kept bit-exact; `decompress_paf` restores the exact bytes. Nothing is written unless `options['output']` is given: then the raw gzip bytes are written to that file (an existing file needs `overwrite=True`).

- **API:** `dnakit.io_formats.compress_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: text (`str`, encoded as UTF-8), `bytes`/`bytearray`, or the path (`str` or `Path`) of an existing file whose bytes are compressed (a directory raises `InputFormatError`; only the zarr functions accept one). A one-line string that is an existing path is read as a path. A `Path`, or a string with a known file suffix such as `.2bit`, that does not exist raises `InputFormatError` (code `IO_FORMATS_FILE_NOT_FOUND`) instead of being compressed as text. Mappings, lists and numbers are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `level`: int 1-9, default 6: gzip compression level
    - `output`: path of a .gz file to write (optional; nothing is written otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `encoding` (`"gzip+base64"`), `data` (the base64 text), `original_bytes` (size before compression), `compressed_bytes` (size of the gzip stream), `sha256_original` (SHA-256 of the uncompressed bytes) and `output_path` when `options['output']` is given. Unusable inputs raise `InputFormatError`; invalid or unknown options (for example `level` outside 1-9), an existing output file without `overwrite=True` and inputs above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import compress_paf

envelope = compress_paf('q1\t100\t0\t100\t+\tt1\t500\t10\t110\t95\t100\t60\n', verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAAys05DQ0MOA0AJPanCWGnKYGIA6nIRBbmoKFzQy4AKXfwIgnAAAA",
    "original_bytes": 39,
    "compressed_bytes": 48,
    "sha256_original": "b5f356a5606c2b16c7d2b21441ca4b368d794e6aa5cb9dc4cd87d84a373fcbc6"
  }
}
```

## 9) compress_ped

- **Purpose:** Compress a PLINK `.ped` genotype table (text, bytes or a file) with deterministic gzip and return it base64-encoded.

- **Method:** The input is turned into bytes (see Input) and compressed with gzip at `options['level']` (1-9, default 6). The gzip header is fixed (modification time 0, no file name), so the same input and level always give the same `data`. The result is base64-encoded and returned with the sizes and the SHA-256 of the uncompressed bytes. The content is never parsed or validated as a PLINK `.ped` genotype table, and binary content is kept bit-exact; `decompress_ped` restores the exact bytes. Nothing is written unless `options['output']` is given: then the raw gzip bytes are written to that file (an existing file needs `overwrite=True`).

- **API:** `dnakit.io_formats.compress_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: text (`str`, encoded as UTF-8), `bytes`/`bytearray`, or the path (`str` or `Path`) of an existing file whose bytes are compressed (a directory raises `InputFormatError`; only the zarr functions accept one). A one-line string that is an existing path is read as a path. A `Path`, or a string with a known file suffix such as `.2bit`, that does not exist raises `InputFormatError` (code `IO_FORMATS_FILE_NOT_FOUND`) instead of being compressed as text. Mappings, lists and numbers are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `level`: int 1-9, default 6: gzip compression level
    - `output`: path of a .gz file to write (optional; nothing is written otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `encoding` (`"gzip+base64"`), `data` (the base64 text), `original_bytes` (size before compression), `compressed_bytes` (size of the gzip stream), `sha256_original` (SHA-256 of the uncompressed bytes) and `output_path` when `options['output']` is given. Unusable inputs raise `InputFormatError`; invalid or unknown options (for example `level` outside 1-9), an existing output file without `overwrite=True` and inputs above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import compress_ped

envelope = compress_ped('F1 I1 0 0 1 -9 A A G G\n', verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAA3MzVPA0VDAAQkMFXUsFRyB0V3DnAgC6+u2xFwAAAA==",
    "original_bytes": 23,
    "compressed_bytes": 40,
    "sha256_original": "e51dd9783e6ea7af725321bed570ce6692918994b6b8d982c24cf918e3d6a858"
  }
}
```

## 10) compress_phylip

- **Purpose:** Compress a PHYLIP alignment (text, bytes or a file) with deterministic gzip and return it base64-encoded.

- **Method:** The input is turned into bytes (see Input) and compressed with gzip at `options['level']` (1-9, default 6). The gzip header is fixed (modification time 0, no file name), so the same input and level always give the same `data`. The result is base64-encoded and returned with the sizes and the SHA-256 of the uncompressed bytes. The content is never parsed or validated as a PHYLIP alignment, and binary content is kept bit-exact; `decompress_phylip` restores the exact bytes. Nothing is written unless `options['output']` is given: then the raw gzip bytes are written to that file (an existing file needs `overwrite=True`).

- **API:** `dnakit.io_formats.compress_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: text (`str`, encoded as UTF-8), `bytes`/`bytearray`, or the path (`str` or `Path`) of an existing file whose bytes are compressed (a directory raises `InputFormatError`; only the zarr functions accept one). A one-line string that is an existing path is read as a path. A `Path`, or a string with a known file suffix such as `.2bit`, that does not exist raises `InputFormatError` (code `IO_FORMATS_FILE_NOT_FOUND`) instead of being compressed as text. Mappings, lists and numbers are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `level`: int 1-9, default 6: gzip compression level
    - `output`: path of a .gz file to write (optional; nothing is written otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `encoding` (`"gzip+base64"`), `data` (the base64 text), `original_bytes` (size before compression), `compressed_bytes` (size of the gzip stream), `sha256_original` (SHA-256 of the uncompressed bytes) and `output_path` when `options['output']` is given. Unusable inputs raise `InputFormatError`; invalid or unknown options (for example `level` outside 1-9), an existing output file without `overwrite=True` and inputs above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import compress_phylip

envelope = compress_phylip(' 2 4\nA ACGT\nB ACGA\n', verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAA1MwUjDhclRwdHYP4XICUY5cAGNU7OsTAAAA",
    "original_bytes": 19,
    "compressed_bytes": 36,
    "sha256_original": "64de1e5a4b3f8ab13596909fb54c08ddde79434041d40df2eb2b4ed435815d7d"
  }
}
```

## 11) compress_twobit

- **Purpose:** Compress a UCSC 2bit sequence file (text, bytes or a file) with deterministic gzip and return it base64-encoded.

- **Method:** The input is turned into bytes (see Input) and compressed with gzip at `options['level']` (1-9, default 6). The gzip header is fixed (modification time 0, no file name), so the same input and level always give the same `data`. The result is base64-encoded and returned with the sizes and the SHA-256 of the uncompressed bytes. The content is never parsed or validated as a UCSC 2bit sequence file, and binary content is kept bit-exact; `decompress_twobit` restores the exact bytes. Nothing is written unless `options['output']` is given: then the raw gzip bytes are written to that file (an existing file needs `overwrite=True`).

- **API:** `dnakit.io_formats.compress_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: text (`str`, encoded as UTF-8), `bytes`/`bytearray`, or the path (`str` or `Path`) of an existing file whose bytes are compressed (a directory raises `InputFormatError`; only the zarr functions accept one). A one-line string that is an existing path is read as a path. A `Path`, or a string with a known file suffix such as `.2bit`, that does not exist raises `InputFormatError` (code `IO_FORMATS_FILE_NOT_FOUND`) instead of being compressed as text. Mappings, lists and numbers are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `level`: int 1-9, default 6: gzip compression level
    - `output`: path of a .gz file to write (optional; nothing is written otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `encoding` (`"gzip+base64"`), `data` (the base64 text), `original_bytes` (size before compression), `compressed_bytes` (size of the gzip stream), `sha256_original` (SHA-256 of the uncompressed bytes) and `output_path` when `options['output']` is given. Unusable inputs raise `InputFormatError`; invalid or unknown options (for example `level` outside 1-9), an existing output file without `overwrite=True` and inputs above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import base64
import json

from dnakit.io_formats import compress_twobit, write_twobit

blob = base64.b64decode(write_twobit({"chr1": "ACGTACGTNN"}, verbose=False)["result"]["data_base64"])
envelope = compress_twobit(blob, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAA3NWd5RiAAJGBghgSc4oMpQEMrigghxAzMSAAHPmMAAAOxL/rzQAAAA=",
    "original_bytes": 52,
    "compressed_bytes": 50,
    "sha256_original": "9aa5d8bef92c1ae93359b6716557ed23609d3e7f4719ad399e765e40cb2475d1"
  }
}
```

## 12) compress_zarr

- **Purpose:** Compress a Zarr store (a directory, a file, text or bytes) with deterministic gzip and return it base64-encoded.

- **Method:** The input is turned into bytes (see Input) and compressed with gzip at `options['level']` (1-9, default 6). The gzip header is fixed (modification time 0, no file name), so the same input and level always give the same `data`. The result is base64-encoded and returned with the sizes and the SHA-256 of the uncompressed bytes. The content is never parsed or validated as a Zarr store, and binary content is kept bit-exact; `decompress_zarr` restores the exact bytes. Nothing is written unless `options['output']` is given: then the raw gzip bytes are written to that file (an existing file needs `overwrite=True`). A directory is archived first: its files and subdirectories (paths relative to the directory) are written to a tar archive in sorted order with times and owners zeroed, so the archive is reproducible, and the tar archive is what gets compressed. Symbolic links and special files inside the directory raise `InputFormatError`.

- **API:** `dnakit.io_formats.compress_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: text (`str`, encoded as UTF-8), `bytes`/`bytearray`, or the path (`str` or `Path`) of an existing file whose bytes are compressed or of a directory (see Method). A one-line string that is an existing path is read as a path. A `Path`, or a string with a known file suffix such as `.2bit`, that does not exist raises `InputFormatError` (code `IO_FORMATS_FILE_NOT_FOUND`) instead of being compressed as text. Mappings, lists and numbers are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `level`: int 1-9, default 6: gzip compression level
    - `output`: path of a .gz file to write (optional; nothing is written otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `encoding` (`"gzip+base64"`), `data` (the base64 text), `original_bytes` (size before compression), `compressed_bytes` (size of the gzip stream), `sha256_original` (SHA-256 of the uncompressed bytes), `archive` (`"tar"`, only when the input was a directory) and `output_path` when `options['output']` is given. Unusable inputs raise `InputFormatError`; invalid or unknown options (for example `level` outside 1-9), an existing output file without `overwrite=True` and inputs above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json
import shutil
from pathlib import Path

from dnakit.io_formats import compress_zarr

store = Path("demo_store.zarr")
(store / "a").mkdir(parents=True, exist_ok=True)
(store / "zarr.json").write_text('{"zarr_format": 3, "node_type": "group"}')
(store / "a" / "zarr.json").write_text('{"zarr_format": 3, "node_type": "array"}')

envelope = compress_zarr(store, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
shutil.rmtree(store)
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "encoding": "gzip+base64",
    "data": "H4sIAAAAAAAAA+3T0QqCMBiG4V3K2HGklnPQzcggC4KcTD2w6N6bHYUUQTATep+Tf2MH++Hjs4mILg2M1o8ZTOeLc5FttJA6/mpC9G1nvZRzfLVENrlY79en1tXR/hhDLfL8ff56mr8xJhMyjbbRkz/P/6rG/MuD82fbqZ3crqSq3b4qu6Gpwl2FVzuo26/3RBzx2/9V/7Wm/7P42P+jd31D/wEAAAAAAAAAAAAAWLI7rWWEdQAoAAA=",
    "original_bytes": 10240,
    "compressed_bytes": 185,
    "sha256_original": "ae3c6d9a98200511fcb14ad7fa368b2b6fee1b53c073c3b76128f516cb2a8421",
    "archive": "tar"
  }
}
```

## 13) convert

- **Purpose:** Convert a file of any supported format (map, ped, paf, newick, nexus, phylip, twobit, zarr; FASTA as a source) to another format.

- **Method:** The format is taken from `options['format']` (`map`, `newick`, `nexus`, `paf`, `ped`, `phylip`, `twobit` or `zarr`; `fasta` only as the source of `convert`) or detected, first from the path (suffixes `.map`, `.ped`, `.paf`, `.nwk`/`.newick`/`.tre`/`.tree`, `.nex`/`.nexus`/`.nxs`, `.phy`/`.phylip`, `.2bit`, `.fa`/`.fasta`/`.fna`, with a trailing `.gz` ignored; a directory is taken as zarr) and otherwise from the first 64 KiB of content (2bit signature, `#NEXUS`, a leading `(` for Newick, a leading `>` for FASTA, a first line with at least 12 tab-separated columns whose columns 2 to 4 are integers for PAF, an `ntax nchar` header for PHYLIP, four columns ending in an integer for MAP, an even number of at least six columns for PED). The call is then executed by the same engine as the format-specific function (for example `convert_ped`), so the options and the `result` are those of that function. If the format cannot be detected `InputFormatError` (`IO_FORMATS_UNKNOWN_FORMAT`) is raised. The target is `options['to']` (required) and must be one of the targets of the detected format. A FASTA input with any operation other than `convert` raises `ConfigurationError`.

- **API:** `dnakit.io_formats.convert(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a file of any supported format, as accepted by its format-specific function (path, text or bytes). `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Besides the keys below, the options of the format-specific function of the detected format are accepted. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `to`: target format name (required)
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `format`: name of the input format (map, newick, nexus, paf, ped, phylip, twobit, zarr); detected from the suffix/content when omitted

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: exactly the `result` of the format-specific function that was used; its `format` field tells which one (for FASTA input to `convert`, `source_format` is `fasta`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import convert

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = convert(PED, options={"to": "tsv"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "tsv",
    "record_count": 3,
    "source_format": "ped",
    "bytes": 101,
    "sha256": "288b2c84cc53c5a994f59ab2bcc18e128999dacfd1d9ffe7bf1629c930e9ea3f",
    "text": "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tV1\tV2\nF1\tI1\t0\t0\t1\t2\tA/A\tC/G\nF1\tI2\t0\t0\t2\t1\tA/G\t0/0\nF1\tI3\tI1\tI2\t1\t-9\tG/G\tC/C\n"
  }
}
```

## 14) convert_map

- **Purpose:** Convert a PLINK `.map` file to another format (`map`, `bed`, `tsv` or `json`).

- **Method:** `map` writes the normalised tab-separated file; `bed` writes `chrom`, `bp-1`, `bp`, `ID` per variant (BED is 0-based half-open, so a 1-based bp becomes the interval [bp-1, bp)); `tsv` adds the header `chrom id cm pos`; `json` writes the list of records. `options['to']` is required.

- **API:** `dnakit.io_formats.convert_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.map` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `to`: target format name (required): map, bed, tsv, json
    - `output`: path to write (optional)
    - `overwrite`: bool, default False

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format` (the target format), `record_count`, `source_format`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import convert_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = convert_map(MAP, options={"to": "bed"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "bed",
    "record_count": 4,
    "source_format": "map",
    "bytes": 46,
    "sha256": "c3b8e07bd896ef7348d8fd8270d0bf63a95fb78bc265fc8cd8488a90e1240e9f",
    "text": "1\t99\t100\trs1\n1\t49\t50\trs2\n2\t9\t10\trs3\nX\t4\t5\trs4\n"
  }
}
```

## 15) convert_newick

- **Purpose:** Convert a Newick file to another format (`newick`, `nexus`, `tsv` or `json`).

- **Method:** `newick` writes the normalised file; `nexus` writes a TAXA block with the distinct tip names and a TREES block with trees named `tree1`, `tree2`, ...; `tsv` writes the node table `tree, node_id, parent_id, name, length, is_tip` (tree index from 0, empty parent for the root); `json` writes the list of records with flat node tables. `options['to']` is required.

- **API:** `dnakit.io_formats.convert_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a Newick tree file (one or more `;`-terminated trees), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `to`: target format name (required): newick, nexus, tsv, json
    - `output`: path to write (optional)
    - `overwrite`: bool, default False

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format` (the target format), `record_count`, `source_format`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import convert_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = convert_newick(TREES, options={"to": "tsv"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "tsv",
    "record_count": 2,
    "source_format": "newick",
    "bytes": 188,
    "sha256": "45d02b5266f4e4b290b362ea416595353d61720d12bf5f014eba3e7b64903f3f",
    "text": "tree\tnode_id\tparent_id\tname\tlength\tis_tip\n0\t0\t\t\t\t0\n0\t1\t0\t\t0.5\t0\n0\t2\t1\tA\t1\t1\n0\t3\t1\tB\t2\t1\n0\t4\t0\t\t0.5\t0\n0\t5\t4\tC\t1\t1\n0\t6\t4\tD\t1\t1\n1\t0\t\t\t\t0\n1\t1\t0\tA\t\t1\n1\t2\t0\tB\t\t1\n1\t3\t0\t\t\t0\n1\t4\t3\tC\t\t1\n1\t5\t3\tD\t\t1\n"
  }
}
```

## 16) convert_nexus

- **Purpose:** Convert a NEXUS file to another format (`nexus`, `phylip`, `fasta`, `twobit`, `newick` or `json`).

- **Method:** `nexus` writes the normalised DATA block and TREES block; `phylip` and `fasta` write the matrix (taxon names containing whitespace cannot be written to PHYLIP and raise `InputFormatError`); `twobit` packs the nucleotide matrix into a 2bit file (non-ACGT letters become N blocks); `newick` writes one tree per line (a file without trees raises `InputFormatError`); `json` writes an object with `taxon_count`, `site_count`, `datatype`, `trees` and `records`. `options['to']` is required. The `phylip`, `fasta` and `twobit` targets hold only the matrix: when the file also has a TREES block, the matrix is converted and a warning that the trees were dropped is returned.

- **API:** `dnakit.io_formats.convert_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a NEXUS file (`#NEXUS`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `to`: target format name (required): nexus, phylip, fasta, twobit, newick, json
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `layout`: 'sequential' (default) or 'interleaved' for PHYLIP/NEXUS output
    - `block_size`: int, default 60: sites per block of an interleaved output
    - `strict_names`: bool, default False: write strict 10-character PHYLIP names (error when longer)
    - `line_width`: int, default 60: FASTA line width (0 = one line per sequence)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format` (the target format), `record_count`, `source_format`, `bytes` and `sha256` of the serialised output, and the content as `text` (`data_base64` and an extra `sequence_count` for the `twobit` target); with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. When trees were dropped (targets `phylip`, `fasta`, `twobit`), `result` also has `warnings` (a list with one message) and the same message is in the envelope's `warnings`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import convert_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""


envelope = convert_nexus(NEXUS, options={"to": "twobit"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [
    "The 1 tree(s) of the TREES block were dropped: the twobit format holds only the matrix."
  ],
  "result": {
    "format": "twobit",
    "record_count": 3,
    "sequence_count": 3,
    "warnings": [
      "The 1 tree(s) of the TREES block were dropped: the twobit format holds only the matrix."
    ],
    "source_format": "nexus",
    "bytes": 98,
    "sha256": "90833cd90314c8282040eb0e07e8cdb33abbfa3834000cbcbd4ba6e1ad8781e5",
    "data_base64": "QydBGgAAAAADAAAAAAAAAAFBJAAAAAFCNgAAAANDIGRIAAAACAAAAAAAAAAAAAAAAAAAAJycCAAAAAAAAAAAAAAAAAAAAJyeCAAAAAEAAAADAAAAAQAAAAAAAAAAAAAAnJw="
  }
}
```

## 17) convert_paf

- **Purpose:** Convert a PAF file to another format (`paf`, `bed`, `tsv` or `json`).

- **Method:** `paf` writes the normalised file; `bed` writes the target interval per alignment (`target`, `target_start`, `target_end`, `query`, `mapq`, `strand`; no header); `tsv` writes a header with the 12 columns plus one column per tag name (sorted) and the tag values (empty when absent); `json` writes the list of records. `options['to']` is required.

- **API:** `dnakit.io_formats.convert_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PAF (Pairwise mApping Format) file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `to`: target format name (required): paf, bed, tsv, json
    - `output`: path to write (optional)
    - `overwrite`: bool, default False

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format` (the target format), `record_count`, `source_format`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import convert_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = convert_paf(PAF, options={"to": "tsv"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "tsv",
    "record_count": 2,
    "source_format": "paf",
    "bytes": 216,
    "sha256": "9fdb7fe2e77ddbbc6c4e27f7eec4b229553c7255da17988ad1e489b3afbf6752",
    "text": "query_name\tquery_length\tquery_start\tquery_end\tstrand\ttarget_name\ttarget_length\ttarget_start\ttarget_end\tmatches\tblock_length\tmapq\tcg\ttp\nq1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\t90M\tP\nq2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\t\t\n"
  }
}
```

## 18) convert_ped

- **Purpose:** Convert a PLINK `.ped` file to another format (`ped`, `tsv` or `json`).

- **Method:** `ped` writes the normalised tab-separated file; `tsv` writes the header `FID IID PAT MAT SEX PHENO V1..Vk` with genotypes as `a1/a2`; `json` writes an object with `variant_count` and `records`. `options['to']` is required.

- **API:** `dnakit.io_formats.convert_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.ped` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `to`: target format name (required): ped, tsv, json
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `max_cells`: int, default 20000000

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format` (the target format), `record_count`, `source_format`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import convert_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = convert_ped(PED, options={"to": "tsv"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "tsv",
    "record_count": 3,
    "source_format": "ped",
    "bytes": 101,
    "sha256": "288b2c84cc53c5a994f59ab2bcc18e128999dacfd1d9ffe7bf1629c930e9ea3f",
    "text": "FID\tIID\tPAT\tMAT\tSEX\tPHENO\tV1\tV2\nF1\tI1\t0\t0\t1\t2\tA/A\tC/G\nF1\tI2\t0\t0\t2\t1\tA/G\t0/0\nF1\tI3\tI1\tI2\t1\t-9\tG/G\tC/C\n"
  }
}
```

## 19) convert_phylip

- **Purpose:** Convert a PHYLIP file to another format (`phylip`, `nexus`, `fasta`, `twobit` or `json`).

- **Method:** `phylip` re-lays out the alignment (`layout`, `block_size`, `strict_names`); `nexus` writes a DATA block; `fasta` writes records wrapped at `line_width` (0 = one line); `twobit` packs the nucleotide sequences into a 2bit file (non-ACGT letters become N blocks); `json` writes the records. `options['to']` is required.

- **API:** `dnakit.io_formats.convert_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PHYLIP alignment file (header `ntax nchar`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `to`: target format name (required): phylip, nexus, fasta, twobit, json
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `datatype`: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing characters (gap, ?, N for nucleotides, X for protein)
    - `input_layout`: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
    - `input_strict_names`: bool, default False: input names occupy the first 10 characters of each row (otherwise the first whitespace-delimited word)
    - `layout`: 'sequential' (default) or 'interleaved' for PHYLIP/NEXUS output
    - `block_size`: int, default 60: sites per block of an interleaved output
    - `strict_names`: bool, default False: write strict 10-character PHYLIP names (error when longer)
    - `line_width`: int, default 60: FASTA line width (0 = one line per sequence)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format` (the target format), `record_count`, `source_format`, `bytes` and `sha256` of the serialised output, and the content as `text` (`data_base64` and an extra `sequence_count` for the `twobit` target); with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import convert_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = convert_phylip(PHYLIP, options={"to": "fasta"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "fasta",
    "record_count": 3,
    "source_format": "phylip",
    "bytes": 36,
    "sha256": "5f65c0036a72aabb466c7dc332dea148b25953356b317de61c897fe465de5dea",
    "text": ">A\nACGTACGT\n>B\nACGTACGA\n>C\nACG-ACGT\n"
  }
}
```

## 20) convert_twobit

- **Purpose:** Convert a 2bit file to another format (`twobit`, `fasta` or `json`).

- **Method:** `fasta` writes the sequences wrapped at `line_width` (0 = one line), keeping soft-masking as lower case unless `soft_mask=False`; `json` writes a list with `name`, `length`, `n_blocks`, `mask_blocks` and `sequence`; `twobit` re-encodes the file (for example to change `byte_order`). `names` limits the conversion to some sequences; decoding more than `max_bases` bases is refused. FASTA to 2bit is `convert(format='fasta', to='twobit')`. `options['to']` is required.

- **API:** `dnakit.io_formats.convert_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a UCSC 2bit genome file, given as a path (`str` or `Path`) or as its raw `bytes` (a `str` that is not an existing path is rejected, because text cannot hold a binary file); both byte orders are accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `to`: target format name (required): twobit, fasta, json
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `line_width`: int, default 60: FASTA line width (0 = one line)
    - `soft_mask`: bool, default True
    - `max_bases`: int, default 50000000
    - `names`: list of sequence names to convert (default all)
    - `byte_order`: 'little' (default) or 'big': byte order of the written 2bit file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format` (the target format), `record_count`, `source_format`, `bytes` and `sha256` of the serialised output, and the content as `text` (`data_base64` for the `twobit` target, which also has `sequence_count`); with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of the content. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import base64
import json

from dnakit.io_formats import convert_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = convert_twobit(TWOBIT, options={"to": "fasta", "line_width": 5}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "fasta",
    "record_count": 2,
    "source_format": "twobit",
    "bytes": 30,
    "sha256": "d7f761cab7931fa25acad2f48d3b4914cca20d16f8ba89ba150358615295ea28",
    "text": ">s1\nACGTN\nNacgt\n>s2\nGGCCA\nATT\n"
  }
}
```

## 21) convert_zarr

- **Purpose:** Convert a Zarr store to another format (`zarr`, `json`, `tsv` or `npz`).

- **Method:** `json` writes metadata and values; `tsv` writes one 1-D or 2-D array (`array` selects it, required when the store has several) tab-separated without header; `npz` writes a numpy archive of all arrays; `zarr` copies the store to `output`, which is required for this target. Converting more than `max_elements` elements is refused. `options['to']` is required.

- **API:** `dnakit.io_formats.convert_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the path (`str` or `Path`) of a Zarr directory store (zarr v2 or v3; needs the optional `zarr>=3` and `numpy` packages); text content, bytes and gzip are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `to`: target format name (required): zarr, json, tsv, npz
    - `output`: path to write (required for the 'zarr' target, optional otherwise)
    - `overwrite`: bool, default False: replace an existing zarr store at output
    - `array`: array path for the 'tsv' target (1-D or 2-D)
    - `max_elements`: int, default 5000000: refuse to materialise more array elements than this

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: for `json`, `tsv` and `npz` the usual `format` (the target), `record_count`, `bytes`, `sha256` and `text` (`data_base64` for `npz`), or `output_path` instead when `options['output']` is given (an existing file needs `overwrite=True`); for `zarr` `format`, `record_count` and `output_path`. If `zarr` (>= 3) or `numpy` cannot be imported the function returns `status: "not_implemented"` with `required_dependency: "pip install 'zarr>=3' numpy"` and a `message` (`zarr and numpy are required: <import error>`, or `zarr <version> is installed but this adapter needs zarr >= 3.`) instead of a `result`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import convert_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = convert_zarr("demo.zarr", options={"to": "tsv", "array": "dosage"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "tsv",
    "record_count": 4,
    "bytes": 16,
    "sha256": "b593e6e0fd7a0d5575370d49623eabd2848d22f9592cd6b8f9e90695c6ec6d57",
    "text": "0\t1\n1\t1\n0\t0\n2\t1\n"
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

- **Purpose:** Decode gzip data (the `data` of `compress_map`, gzip bytes or a `.gz` file) back into the exact original bytes of the PLINK `.map` marker table.

- **Method:** The input is turned into a gzip stream: bytes starting with the gzip magic (`1f 8b`) are used as they are, anything else must be valid base64 text (whitespace is ignored) of a gzip stream. The stream is decompressed with a limit of `max_bytes` (concatenated gzip members are accepted; truncated or corrupt data and trailing garbage are errors). If the bytes are valid UTF-8 they are returned as `text`, otherwise as `data_base64`, so the round trip is bit-exact; with `options['output']` the bytes are written to that file instead (an existing file needs `overwrite=True`). No PLINK `.map` marker table parsing or validation is done.

- **API:** `dnakit.io_formats.decompress_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the base64 text, i.e. `result["data"]` of `compress_map`; gzip `bytes`; or the path (`str` or `Path`) of a `.gz` file (raw gzip, or a file holding the base64 text). The whole envelope or its `result` mapping is not accepted. A `Path`, or a string with a known file suffix such as `.gz`, that does not exist raises `InputFormatError`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs and decompressed data larger than this
    - `output`: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `original_bytes` (decompressed size), `compressed_bytes`, `sha256` (SHA-256 of the decompressed bytes) and the content as `text` (valid UTF-8) or `data_base64`, or `output_path` when `options['output']` is given. An input that is not valid base64, not gzip, truncated or corrupt raises `InputFormatError` (codes `IO_FORMATS_BAD_BASE64`, `IO_FORMATS_BAD_GZIP`) with a hint; it no longer returns `status: "ok"` with a failure result. Invalid or unknown options, an existing output without `overwrite=True` and data above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import decompress_map

payload = 'H4sIAAAAAAAAAzPkLCo25DTgNDQwMOAyBHKMgBwjEAcAAo9oVBoAAAA='
envelope = decompress_map(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "original_bytes": 26,
    "compressed_bytes": 41,
    "sha256": "0d63a7da7a1761ff406f984298d593a88a953100925de4daf0f0b345438301e9",
    "text": "1\trs1\t0\t1000\n1\trs2\t0\t2000\n"
  }
}
```

## 25) decompress_newick

- **Purpose:** Decode gzip data (the `data` of `compress_newick`, gzip bytes or a `.gz` file) back into the exact original bytes of the Newick tree file.

- **Method:** The input is turned into a gzip stream: bytes starting with the gzip magic (`1f 8b`) are used as they are, anything else must be valid base64 text (whitespace is ignored) of a gzip stream. The stream is decompressed with a limit of `max_bytes` (concatenated gzip members are accepted; truncated or corrupt data and trailing garbage are errors). If the bytes are valid UTF-8 they are returned as `text`, otherwise as `data_base64`, so the round trip is bit-exact; with `options['output']` the bytes are written to that file instead (an existing file needs `overwrite=True`). No Newick tree file parsing or validation is done.

- **API:** `dnakit.io_formats.decompress_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the base64 text, i.e. `result["data"]` of `compress_newick`; gzip `bytes`; or the path (`str` or `Path`) of a `.gz` file (raw gzip, or a file holding the base64 text). The whole envelope or its `result` mapping is not accepted. A `Path`, or a string with a known file suffix such as `.gz`, that does not exist raises `InputFormatError`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs and decompressed data larger than this
    - `output`: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `original_bytes` (decompressed size), `compressed_bytes`, `sha256` (SHA-256 of the decompressed bytes) and the content as `text` (valid UTF-8) or `data_base64`, or `output_path` when `options['output']` is given. An input that is not valid base64, not gzip, truncated or corrupt raises `InputFormatError` (codes `IO_FORMATS_BAD_BASE64`, `IO_FORMATS_BAD_GZIP`) with a hint; it no longer returns `status: "ok"` with a failure result. Invalid or unknown options, an existing output without `overwrite=True` and data above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import decompress_newick

payload = 'H4sIAAAAAAAAA9PQcLQy1HGyMtIEUs5WxprWADsnB7wSAAAA'
envelope = decompress_newick(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "original_bytes": 18,
    "compressed_bytes": 36,
    "sha256": "1687864f64fa679bf85d6ff2ba780918d7c19e3965bcf14f658c9ee9c568c968",
    "text": "((A:1,B:2):1,C:3);"
  }
}
```

## 26) decompress_nexus

- **Purpose:** Decode gzip data (the `data` of `compress_nexus`, gzip bytes or a `.gz` file) back into the exact original bytes of the NEXUS file.

- **Method:** The input is turned into a gzip stream: bytes starting with the gzip magic (`1f 8b`) are used as they are, anything else must be valid base64 text (whitespace is ignored) of a gzip stream. The stream is decompressed with a limit of `max_bytes` (concatenated gzip members are accepted; truncated or corrupt data and trailing garbage are errors). If the bytes are valid UTF-8 they are returned as `text`, otherwise as `data_base64`, so the round trip is bit-exact; with `options['output']` the bytes are written to that file instead (an existing file needs `overwrite=True`). No NEXUS file parsing or validation is done.

- **API:** `dnakit.io_formats.decompress_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the base64 text, i.e. `result["data"]` of `compress_nexus`; gzip `bytes`; or the path (`str` or `Path`) of a `.gz` file (raw gzip, or a file holding the base64 text). The whole envelope or its `result` mapping is not accepted. A `Path`, or a string with a known file suffix such as `.gz`, that does not exist raises `InputFormatError`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs and decompressed data larger than this
    - `output`: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `original_bytes` (decompressed size), `compressed_bytes`, `sha256` (SHA-256 of the decompressed bytes) and the content as `text` (valid UTF-8) or `data_base64`, or `output_path` when `options['output']` is given. An input that is not valid base64, not gzip, truncated or corrupt raises `InputFormatError` (codes `IO_FORMATS_BAD_BASE64`, `IO_FORMATS_BAD_GZIP`) with a hint; it no longer returns `status: "ok"` with a failure result. Invalid or unknown options, an existing output without `overwrite=True` and data above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import decompress_nexus

payload = 'H4sIAAAAAAAAA1P2c40IDeZycnX39FMIcYxwtOZSUHDx9HX1C/b09wtW8AOK2RqBBIEMH0cnV59gBUcFJ2suVz8Xay4AkYWK6z4AAAA='
envelope = decompress_nexus(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "original_bytes": 62,
    "compressed_bytes": 77,
    "sha256": "80a4f6422311b951751d59815bb5a36465fe705e84fb840ecb8f30d8c4a36a8e",
    "text": "#NEXUS\nBEGIN TAXA;\n  DIMENSIONS NTAX=2;\n  TAXLABELS A B;\nEND;\n"
  }
}
```

## 27) decompress_paf

- **Purpose:** Decode gzip data (the `data` of `compress_paf`, gzip bytes or a `.gz` file) back into the exact original bytes of the PAF alignment table.

- **Method:** The input is turned into a gzip stream: bytes starting with the gzip magic (`1f 8b`) are used as they are, anything else must be valid base64 text (whitespace is ignored) of a gzip stream. The stream is decompressed with a limit of `max_bytes` (concatenated gzip members are accepted; truncated or corrupt data and trailing garbage are errors). If the bytes are valid UTF-8 they are returned as `text`, otherwise as `data_base64`, so the round trip is bit-exact; with `options['output']` the bytes are written to that file instead (an existing file needs `overwrite=True`). No PAF alignment table parsing or validation is done.

- **API:** `dnakit.io_formats.decompress_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the base64 text, i.e. `result["data"]` of `compress_paf`; gzip `bytes`; or the path (`str` or `Path`) of a `.gz` file (raw gzip, or a file holding the base64 text). The whole envelope or its `result` mapping is not accepted. A `Path`, or a string with a known file suffix such as `.gz`, that does not exist raises `InputFormatError`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs and decompressed data larger than this
    - `output`: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `original_bytes` (decompressed size), `compressed_bytes`, `sha256` (SHA-256 of the decompressed bytes) and the content as `text` (valid UTF-8) or `data_base64`, or `output_path` when `options['output']` is given. An input that is not valid base64, not gzip, truncated or corrupt raises `InputFormatError` (codes `IO_FORMATS_BAD_BASE64`, `IO_FORMATS_BAD_GZIP`) with a hint; it no longer returns `status: "ok"` with a failure result. Invalid or unknown options, an existing output without `overwrite=True` and data above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import decompress_paf

payload = 'H4sIAAAAAAAAAys05DQ0MOA0AJPanCWGnKYGIA6nIRBbmoKFzQy4AKXfwIgnAAAA'
envelope = decompress_paf(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "original_bytes": 39,
    "compressed_bytes": 48,
    "sha256": "b5f356a5606c2b16c7d2b21441ca4b368d794e6aa5cb9dc4cd87d84a373fcbc6",
    "text": "q1\t100\t0\t100\t+\tt1\t500\t10\t110\t95\t100\t60\n"
  }
}
```

## 28) decompress_ped

- **Purpose:** Decode gzip data (the `data` of `compress_ped`, gzip bytes or a `.gz` file) back into the exact original bytes of the PLINK `.ped` genotype table.

- **Method:** The input is turned into a gzip stream: bytes starting with the gzip magic (`1f 8b`) are used as they are, anything else must be valid base64 text (whitespace is ignored) of a gzip stream. The stream is decompressed with a limit of `max_bytes` (concatenated gzip members are accepted; truncated or corrupt data and trailing garbage are errors). If the bytes are valid UTF-8 they are returned as `text`, otherwise as `data_base64`, so the round trip is bit-exact; with `options['output']` the bytes are written to that file instead (an existing file needs `overwrite=True`). No PLINK `.ped` genotype table parsing or validation is done.

- **API:** `dnakit.io_formats.decompress_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the base64 text, i.e. `result["data"]` of `compress_ped`; gzip `bytes`; or the path (`str` or `Path`) of a `.gz` file (raw gzip, or a file holding the base64 text). The whole envelope or its `result` mapping is not accepted. A `Path`, or a string with a known file suffix such as `.gz`, that does not exist raises `InputFormatError`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs and decompressed data larger than this
    - `output`: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `original_bytes` (decompressed size), `compressed_bytes`, `sha256` (SHA-256 of the decompressed bytes) and the content as `text` (valid UTF-8) or `data_base64`, or `output_path` when `options['output']` is given. An input that is not valid base64, not gzip, truncated or corrupt raises `InputFormatError` (codes `IO_FORMATS_BAD_BASE64`, `IO_FORMATS_BAD_GZIP`) with a hint; it no longer returns `status: "ok"` with a failure result. Invalid or unknown options, an existing output without `overwrite=True` and data above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import decompress_ped

payload = 'H4sIAAAAAAAAA3MzVPA0VDAAQkMFXUsFRyB0V3DnAgC6+u2xFwAAAA=='
envelope = decompress_ped(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "original_bytes": 23,
    "compressed_bytes": 40,
    "sha256": "e51dd9783e6ea7af725321bed570ce6692918994b6b8d982c24cf918e3d6a858",
    "text": "F1 I1 0 0 1 -9 A A G G\n"
  }
}
```

## 29) decompress_phylip

- **Purpose:** Decode gzip data (the `data` of `compress_phylip`, gzip bytes or a `.gz` file) back into the exact original bytes of the PHYLIP alignment.

- **Method:** The input is turned into a gzip stream: bytes starting with the gzip magic (`1f 8b`) are used as they are, anything else must be valid base64 text (whitespace is ignored) of a gzip stream. The stream is decompressed with a limit of `max_bytes` (concatenated gzip members are accepted; truncated or corrupt data and trailing garbage are errors). If the bytes are valid UTF-8 they are returned as `text`, otherwise as `data_base64`, so the round trip is bit-exact; with `options['output']` the bytes are written to that file instead (an existing file needs `overwrite=True`). No PHYLIP alignment parsing or validation is done.

- **API:** `dnakit.io_formats.decompress_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the base64 text, i.e. `result["data"]` of `compress_phylip`; gzip `bytes`; or the path (`str` or `Path`) of a `.gz` file (raw gzip, or a file holding the base64 text). The whole envelope or its `result` mapping is not accepted. A `Path`, or a string with a known file suffix such as `.gz`, that does not exist raises `InputFormatError`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs and decompressed data larger than this
    - `output`: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `original_bytes` (decompressed size), `compressed_bytes`, `sha256` (SHA-256 of the decompressed bytes) and the content as `text` (valid UTF-8) or `data_base64`, or `output_path` when `options['output']` is given. An input that is not valid base64, not gzip, truncated or corrupt raises `InputFormatError` (codes `IO_FORMATS_BAD_BASE64`, `IO_FORMATS_BAD_GZIP`) with a hint; it no longer returns `status: "ok"` with a failure result. Invalid or unknown options, an existing output without `overwrite=True` and data above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import decompress_phylip

payload = 'H4sIAAAAAAAAA1MwUjDhclRwdHYP4XICUY5cAGNU7OsTAAAA'
envelope = decompress_phylip(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "original_bytes": 19,
    "compressed_bytes": 36,
    "sha256": "64de1e5a4b3f8ab13596909fb54c08ddde79434041d40df2eb2b4ed435815d7d",
    "text": " 2 4\nA ACGT\nB ACGA\n"
  }
}
```

## 30) decompress_twobit

- **Purpose:** Decode gzip data (the `data` of `compress_twobit`, gzip bytes or a `.gz` file) back into the exact original bytes of the UCSC 2bit sequence file.

- **Method:** The input is turned into a gzip stream: bytes starting with the gzip magic (`1f 8b`) are used as they are, anything else must be valid base64 text (whitespace is ignored) of a gzip stream. The stream is decompressed with a limit of `max_bytes` (concatenated gzip members are accepted; truncated or corrupt data and trailing garbage are errors). If the bytes are valid UTF-8 they are returned as `text`, otherwise as `data_base64`, so the round trip is bit-exact; with `options['output']` the bytes are written to that file instead (an existing file needs `overwrite=True`). No UCSC 2bit sequence file parsing or validation is done.

- **API:** `dnakit.io_formats.decompress_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the base64 text, i.e. `result["data"]` of `compress_twobit`; gzip `bytes`; or the path (`str` or `Path`) of a `.gz` file (raw gzip, or a file holding the base64 text). The whole envelope or its `result` mapping is not accepted. A `Path`, or a string with a known file suffix such as `.gz`, that does not exist raises `InputFormatError`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs and decompressed data larger than this
    - `output`: path of a file to write the decompressed bytes to (optional; they are returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `original_bytes` (decompressed size), `compressed_bytes`, `sha256` (SHA-256 of the decompressed bytes) and the content as `text` (valid UTF-8) or `data_base64`, or `output_path` when `options['output']` is given. An input that is not valid base64, not gzip, truncated or corrupt raises `InputFormatError` (codes `IO_FORMATS_BAD_BASE64`, `IO_FORMATS_BAD_GZIP`) with a hint; it no longer returns `status: "ok"` with a failure result. Invalid or unknown options, an existing output without `overwrite=True` and data above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import decompress_twobit

payload = 'H4sIAAAAAAAAA3NWd5RiAAJGBghgSc4oMpQEMrigghxAzMSAAHPmMAAAOxL/rzQAAAA='
envelope = decompress_twobit(payload, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "original_bytes": 52,
    "compressed_bytes": 50,
    "sha256": "9aa5d8bef92c1ae93359b6716557ed23609d3e7f4719ad399e765e40cb2475d1",
    "data_base64": "QydBGgAAAAABAAAAAAAAAARjaHIxGQAAAAoAAAABAAAACAAAAAIAAAAAAAAAAAAAAJycAA=="
  }
}
```

## 31) decompress_zarr

- **Purpose:** Decode gzip data (the `data` of `compress_zarr`, gzip bytes or a `.gz` file) back into the exact original bytes of the Zarr store (extracting a directory archive to `options['output']`).

- **Method:** The input is turned into a gzip stream: bytes starting with the gzip magic (`1f 8b`) are used as they are, anything else must be valid base64 text (whitespace is ignored) of a gzip stream. The stream is decompressed with a limit of `max_bytes` (concatenated gzip members are accepted; truncated or corrupt data and trailing garbage are errors). If the bytes are valid UTF-8 they are returned as `text`, otherwise as `data_base64`, so the round trip is bit-exact; with `options['output']` the bytes are written to that file instead (an existing file needs `overwrite=True`). No Zarr store parsing or validation is done. If the decompressed bytes are a tar archive (what `compress_zarr` makes from a directory), `options['output']` names the directory to extract it to (it must be new or empty unless `overwrite=True`) and the archive is extracted safely: all members are checked before anything is written, and absolute paths, `..` components, links, devices and other special members are refused with `InputFormatError` (code `IO_FORMATS_UNSAFE_ARCHIVE`). Without `output` nothing is written and the tar archive is returned as `data_base64`.

- **API:** `dnakit.io_formats.decompress_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the base64 text, i.e. `result["data"]` of `compress_zarr`; gzip `bytes`; or the path (`str` or `Path`) of a `.gz` file (raw gzip, or a file holding the base64 text). The whole envelope or its `result` mapping is not accepted. A `Path`, or a string with a known file suffix such as `.gz`, that does not exist raises `InputFormatError`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs and decompressed data larger than this
    - `output`: a directory to extract an archive into (it must be empty or new unless overwrite=True), or the path of a file for any other data (optional; without it nothing is written)
    - `overwrite`: bool, default False: allow a non-empty directory (existing files are replaced) or an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `original_bytes` (decompressed size), `compressed_bytes`, `sha256` (SHA-256 of the decompressed bytes) and the content as `text` (valid UTF-8) or `data_base64`, or `output_path` when `options['output']` is given. For a tar archive `result` has `archive: "tar"` and either `output_path` with `file_count` (extracted files) or `data_base64` (the archive). An input that is not valid base64, not gzip, truncated or corrupt raises `InputFormatError` (codes `IO_FORMATS_BAD_BASE64`, `IO_FORMATS_BAD_GZIP`) with a hint; it no longer returns `status: "ok"` with a failure result. Invalid or unknown options, an existing output without `overwrite=True` and data above `max_bytes` raise `ConfigurationError`.

- **Example code:**

```python
import json
import shutil

from dnakit.io_formats import decompress_zarr

payload = 'H4sIAAAAAAAAA+3T0QqCMBiG4V3K2HGklnPQzcggC4KcTD2w6N6bHYUUQTATep+Tf2MH++Hjs4mILg2M1o8ZTOeLc5FttJA6/mpC9G1nvZRzfLVENrlY79en1tXR/hhDLfL8ff56mr8xJhMyjbbRkz/P/6rG/MuD82fbqZ3crqSq3b4qu6Gpwl2FVzuo26/3RBzx2/9V/7Wm/7P42P+jd31D/wEAAAAAAAAAAAAAWLI7rWWEdQAoAAA='
envelope = decompress_zarr(payload, options={"output": "restored_store"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
shutil.rmtree("restored_store")
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "original_bytes": 10240,
    "compressed_bytes": 185,
    "sha256": "ae3c6d9a98200511fcb14ad7fa368b2b6fee1b53c073c3b76128f516cb2a8421",
    "archive": "tar",
    "file_count": 2,
    "output_path": "restored_store"
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

- **Purpose:** Filter a file of any supported format with the criteria of its format-specific filter function.

- **Method:** The format is taken from `options['format']` (`map`, `newick`, `nexus`, `paf`, `ped`, `phylip`, `twobit` or `zarr`; `fasta` only as the source of `convert`) or detected, first from the path (suffixes `.map`, `.ped`, `.paf`, `.nwk`/`.newick`/`.tre`/`.tree`, `.nex`/`.nexus`/`.nxs`, `.phy`/`.phylip`, `.2bit`, `.fa`/`.fasta`/`.fna`, with a trailing `.gz` ignored; a directory is taken as zarr) and otherwise from the first 64 KiB of content (2bit signature, `#NEXUS`, a leading `(` for Newick, a leading `>` for FASTA, a first line with at least 12 tab-separated columns whose columns 2 to 4 are integers for PAF, an `ntax nchar` header for PHYLIP, four columns ending in an integer for MAP, an even number of at least six columns for PED). The call is then executed by the same engine as the format-specific function (for example `filter_paf`), so the options and the `result` are those of that function. If the format cannot be detected `InputFormatError` (`IO_FORMATS_UNKNOWN_FORMAT`) is raised. A FASTA input with any operation other than `convert` raises `ConfigurationError`.

- **API:** `dnakit.io_formats.filter(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a file of any supported format, as accepted by its format-specific function (path, text or bytes). `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Besides the keys below, the options of the format-specific function of the detected format are accepted. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `format`: name of the input format (map, newick, nexus, paf, ped, phylip, twobit, zarr); detected from the suffix/content when omitted

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: exactly the `result` of the format-specific function that was used; its `format` field tells which one (for FASTA input to `convert`, `source_format` is `fasta`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import filter

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = filter(PAF, options={"min_mapq": 30}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "record_count": 1,
    "input_count": 2,
    "output_count": 1,
    "removed_count": 1,
    "bytes": 55,
    "sha256": "5f960e3807baf139ee66a38539d822ab8784cde4fe8a933cdbb7c1d941d2231b",
    "text": "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
  }
}
```

## 37) filter_map

- **Purpose:** Filter a PLINK `.map` file by chromosome, variant ID, position and genetic distance, and return the filtered `.map` content.

- **Method:** A variant is kept when it passes every given criterion (logical AND): chromosome in/not in the given lists, ID in/not in the given lists, `region` (`chrom:start-end`, 1-based inclusive, same chromosome and `start <= bp <= end`), `min_position`/`max_position` on bp and `min_cm`/`max_cm` on the genetic distance. The file order is kept; without any criterion all variants pass.

- **API:** `dnakit.io_formats.filter_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.map` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `chromosomes`: list of chromosome names to keep
    - `exclude_chromosomes`: list of chromosome names to drop
    - `ids`: list of variant IDs to keep
    - `exclude_ids`: list of variant IDs to drop
    - `region`: 'chrom:start-end' keep variants with start <= bp <= end (1-based, inclusive)
    - `min_position`: int: keep bp >= value
    - `max_position`: int: keep bp <= value
    - `min_cm`: float: keep genetic distance >= value
    - `max_cm`: float: keep genetic distance <= value

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `removed_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import filter_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = filter_map(MAP, options={"chromosomes": ["1", "X"], "max_position": 90}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "record_count": 2,
    "input_count": 4,
    "output_count": 2,
    "removed_count": 2,
    "bytes": 23,
    "sha256": "7175112a87148a8e643383f0e65de0b0832f3778411172821c9bf7c5e8da8881",
    "text": "1\trs2\t0.5\t50\nX\trs4\t0\t5\n"
  }
}
```

## 38) filter_newick

- **Purpose:** Filter a Newick file by tip names, short branches and tip counts, and return the filtered Newick content.

- **Method:** Per tree, in this order: `keep_tips`/`drop_tips` prune the tree to the wanted tips (a tree with no remaining tip is dropped; internal nodes left with one child are merged and their branch lengths added); `collapse_below` contracts every internal branch shorter than the threshold into a polytomy; then `min_tips`/`max_tips` drop trees whose remaining tip count is outside the range.

- **API:** `dnakit.io_formats.filter_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a Newick tree file (one or more `;`-terminated trees), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `keep_tips`: list of tip names to keep (other tips are pruned; unary nodes are merged, summing branch lengths)
    - `drop_tips`: list of tip names to prune
    - `collapse_below`: float: contract internal branches shorter than this into polytomies
    - `min_tips`: int: drop trees with fewer tips (evaluated after pruning)
    - `max_tips`: int: drop trees with more tips (evaluated after pruning)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `removed_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import filter_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = filter_newick(TREES, options={"drop_tips": ["D"]}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "record_count": 2,
    "input_count": 2,
    "output_count": 2,
    "removed_count": 0,
    "bytes": 32,
    "sha256": "f2260e09bd860329f784c631272292e6da494b022ae3ee88587218dff8d2e942",
    "text": "((A:1,B:2):0.5,C:1.5);\n(A,B,C);\n"
  }
}
```

## 39) filter_nexus

- **Purpose:** Filter a NEXUS file by taxon names and per-taxon / per-site missing data, pruning the trees to the kept taxa.

- **Method:** Taxa are kept when they pass every given criterion: in `keep`, not in `exclude`, matching `name_pattern` (`re.search`, at most 256 characters) and with a gap/missing/ambiguous fraction not above `max_missing_fraction`. Then alignment columns whose missing fraction over the kept taxa exceeds `max_site_missing_fraction` are removed. Taxon filters apply to the matrix, the TAXA labels and, by pruning tips, every tree (a tree with no remaining tip is dropped).

- **API:** `dnakit.io_formats.filter_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a NEXUS file (`#NEXUS`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `keep`: list of taxon names to keep
    - `exclude`: list of taxon names to drop
    - `name_pattern`: regular expression (re.search, max 256 characters) a taxon name must match
    - `max_missing_fraction`: float in [0,1]: drop taxa whose gap/missing/ambiguous fraction exceeds it
    - `max_site_missing_fraction`: float in [0,1]: drop alignment columns whose missing fraction (over the kept taxa) exceeds it

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `removed_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import filter_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = filter_nexus(NEXUS, options={"exclude": ["B"]}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "record_count": 3,
    "input_count": 4,
    "output_count": 3,
    "removed_count": 1,
    "bytes": 207,
    "sha256": "a4c11f842c4c2c54c81268d8380d449f4dfe72e33af3929d45d73aa6c735684e",
    "text": "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=2 NCHAR=8;\n    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n    MATRIX\n    A      ACGTACGT\n    'C d'  ACG-ACGT\n    ;\nEND;\n\nBEGIN TREES;\n    TREE t1 = [&R] (A:2,'C d':2);\nEND;\n"
  }
}
```

## 40) filter_paf

- **Purpose:** Filter a PAF file by mapping quality, identity, length, strand, name, region and tag criteria, and return the filtered PAF content.

- **Method:** An alignment is kept when it passes every given criterion (logical AND): `min_mapq`, `min_matches`, `min_block_length`, `min_identity` (matches / block_length), `min_query_coverage` ((query_end - query_start) / query_length), `strand`, query/target name lists to keep or drop, `target_region`/`query_region` (`name:start-end`, kept when the alignment overlaps the 0-based half-open interval on that sequence) and `tags` (every listed tag must be present with that value, compared as text). The input order is kept.

- **API:** `dnakit.io_formats.filter_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PAF (Pairwise mApping Format) file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `min_mapq`: int: keep mapping quality >= value
    - `min_identity`: float: keep matches/block_length >= value
    - `min_block_length`: int: keep alignment block length >= value
    - `min_matches`: int: keep matching bases >= value
    - `min_query_coverage`: float: keep (query_end-query_start)/query_length >= value
    - `strand`: '+' or '-'
    - `query_names`: list of query names to keep
    - `exclude_query_names`: list to drop
    - `target_names`: list of target names to keep
    - `exclude_target_names`: list to drop
    - `target_region`: 'name:start-end' keep alignments overlapping the 0-based half-open target interval
    - `query_region`: 'name:start-end' keep alignments overlapping the 0-based half-open query interval
    - `tags`: mapping {tag: value}: keep alignments carrying every tag with that value (e.g. {'tp': 'P'})

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `removed_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import filter_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = filter_paf(PAF, options={"min_mapq": 30, "tags": {"tp": "P"}}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "record_count": 1,
    "input_count": 2,
    "output_count": 1,
    "removed_count": 1,
    "bytes": 55,
    "sha256": "5f960e3807baf139ee66a38539d822ab8784cde4fe8a933cdbb7c1d941d2231b",
    "text": "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
  }
}
```

## 41) filter_ped

- **Purpose:** Filter a PLINK `.ped` file by family, individual, sex, phenotype and missing-genotype criteria, and return the filtered `.ped` content.

- **Method:** Individual filters are combined with AND: families/individuals to keep or drop, `sex` (1/`male`, 2/`female`, 0/`unknown`), `phenotype` class (`case` = 2, `control` = 1, `missing` = 0 or -9, `observed` = not missing, `quantitative` = any other value), `founders_only` (father and mother both `0`) and `max_missing_rate` (fraction of missing genotype calls per individual). Variant filters then drop genotype columns from every kept row: `variant_indices` keeps the listed 0-based columns, and `max_variant_missing_rate` drops columns whose missing-call fraction, computed over the individuals that passed the other filters, exceeds the threshold. The same column selection has to be applied to the matching `.map` file by the caller.

- **API:** `dnakit.io_formats.filter_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.ped` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `families`: list of family IDs to keep
    - `exclude_families`: list of family IDs to drop
    - `individuals`: list of individual IDs to keep
    - `exclude_individuals`: list of individual IDs to drop
    - `sex`: keep one sex: 1/'male', 2/'female' or 0/'unknown'
    - `phenotype`: keep one class: 'case' (2), 'control' (1), 'missing' (0 or -9), 'observed' (not missing) or 'quantitative' (any other value)
    - `founders_only`: bool, default False: keep individuals whose father and mother IDs are both 0
    - `max_missing_rate`: float in [0,1]: drop individuals whose fraction of missing genotypes exceeds it
    - `variant_indices`: list of 0-based variant columns to keep (apply the same to the .map file)
    - `max_variant_missing_rate`: float in [0,1]: drop variant columns whose missing-call fraction exceeds it (computed over the individuals that pass the other filters)
    - `max_cells`: int, default 20000000

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `removed_count`, `variant_count`, `kept_variant_indices` (the 0-based variant columns that were kept), `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import filter_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = filter_ped(PED, options={"founders_only": True, "max_variant_missing_rate": 0.4}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "record_count": 2,
    "variant_count": 1,
    "input_count": 3,
    "output_count": 2,
    "removed_count": 1,
    "kept_variant_indices": [
      0
    ],
    "bytes": 36,
    "sha256": "8c55765c8c50404f7b360d6c86a60f6191faa94deb8837c14f9e6b013f6d0f56",
    "text": "F1\tI1\t0\t0\t1\t2\tA\tA\nF1\tI2\t0\t0\t2\t1\tA\tG\n"
  }
}
```

## 42) filter_phylip

- **Purpose:** Filter a PHYLIP file by taxon names and per-taxon / per-site missing data.

- **Method:** Taxa are kept when they pass every given criterion: in `keep`, not in `exclude`, matching `name_pattern` (`re.search`, at most 256 characters) and with a gap/missing/ambiguous fraction not above `max_missing_fraction`. Then alignment columns whose missing fraction over the kept taxa exceeds `max_site_missing_fraction` are removed.

- **API:** `dnakit.io_formats.filter_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PHYLIP alignment file (header `ntax nchar`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `datatype`: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing characters (gap, ?, N for nucleotides, X for protein)
    - `input_layout`: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
    - `input_strict_names`: bool, default False: input names occupy the first 10 characters of each row (otherwise the first whitespace-delimited word)
    - `keep`: list of taxon names to keep
    - `exclude`: list of taxon names to drop
    - `name_pattern`: regular expression (re.search, max 256 characters) a taxon name must match
    - `max_missing_fraction`: float in [0,1]: drop taxa whose gap/missing/ambiguous fraction exceeds it
    - `max_site_missing_fraction`: float in [0,1]: drop alignment columns whose missing fraction (over the kept taxa) exceeds it

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `removed_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import filter_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = filter_phylip(PHYLIP, options={"max_site_missing_fraction": 0.2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "record_count": 3,
    "input_count": 3,
    "output_count": 3,
    "removed_count": 0,
    "bytes": 62,
    "sha256": "89eeb0f0715875534c807a2240cb826134cee6d673b22d4197eb4d981a042300",
    "text": " 3 7\nA          ACGACGT\nB          ACGACGA\nC          ACGACGT\n"
  }
}
```

## 43) filter_twobit

- **Purpose:** Filter a 2bit file by sequence name, length and N fraction, and return the filtered 2bit content.

- **Method:** Whole sequences are kept (their records are copied bit-exactly, without re-encoding) when they pass every given criterion: in `names`, not in `exclude`, matching `name_pattern` (`re.search`, at most 256 characters), length between `min_length` and `max_length`, and N fraction (N bases / length) not above `max_n_fraction`.

- **API:** `dnakit.io_formats.filter_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a UCSC 2bit genome file, given as a path (`str` or `Path`) or as its raw `bytes` (a `str` that is not an existing path is rejected, because text cannot hold a binary file); both byte orders are accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `names`: list of sequence names to keep
    - `exclude`: list of names to drop
    - `name_pattern`: regular expression (re.search, max 256 characters) the name must match
    - `min_length`: int: keep sequences with at least this many bases
    - `max_length`: int: at most this many bases
    - `max_n_fraction`: float in [0,1]: drop sequences whose N fraction exceeds it
    - `byte_order`: 'little' (default) or 'big': byte order of the written 2bit file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `removed_count`, `sequence_count`, `bytes` and `sha256` of the serialised output, and the content as `data_base64`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `data_base64`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import base64
import json

from dnakit.io_formats import filter_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = filter_twobit(TWOBIT, options={"min_length": 9}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "record_count": 1,
    "sequence_count": 1,
    "input_count": 2,
    "output_count": 1,
    "removed_count": 1,
    "bytes": 58,
    "sha256": "2c8722ce953b5a3db475f68d7f8f16d9f09381af10b7e1472ac129e3eed9f6cf",
    "data_base64": "QydBGgAAAAABAAAAAAAAAAJzMRcAAAAKAAAAAQAAAAQAAAACAAAAAQAAAAYAAAAEAAAAAAAAAJwJwA=="
  }
}
```

## 44) filter_zarr

- **Purpose:** Filter a Zarr store by array paths and rows along the first axis, writing a new store.

- **Method:** Copies the arrays selected by `arrays` (paths or group prefixes; default all) minus `exclude` into a new store at `output` (required), and with `rows` only the chosen rows along axis 0 (`'start:stop[:step]'`, `{start, stop, step}` or a list of indices; all selected arrays with at least one dimension must have the same length on axis 0). Copying more than `max_elements` elements is refused.

- **API:** `dnakit.io_formats.filter_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the path (`str` or `Path`) of a Zarr directory store (zarr v2 or v3; needs the optional `zarr>=3` and `numpy` packages); text content, bytes and gzip are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `output`: directory path of the new zarr store (required)
    - `overwrite`: bool, default False: replace an existing zarr store at output
    - `arrays`: array paths (or group prefixes) to keep
    - `exclude`: array paths (or group prefixes) to drop
    - `rows`: rows to keep along axis 0: 'start:stop[:step]', {start, stop, step} or a list of indices (all selected arrays with at least one dimension must have the same length on axis 0)
    - `max_elements`: int, default 5000000: refuse to materialise more array elements than this

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count` (number of arrays written), `kept_rows` (when `rows` is given) and `output_path` (no content is returned). If `zarr` (>= 3) or `numpy` cannot be imported the function returns `status: "not_implemented"` with `required_dependency: "pip install 'zarr>=3' numpy"` and a `message` (`zarr and numpy are required: <import error>`, or `zarr <version> is installed but this adapter needs zarr >= 3.`) instead of a `result`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import filter_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = filter_zarr("demo.zarr", options={"output": "filtered.zarr", "overwrite": True, "arrays": ["dosage"], "rows": "1:3"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "record_count": 1,
    "kept_rows": 2,
    "output_path": "filtered.zarr"
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

- **Purpose:** Build a lookup index of a file of any supported format.

- **Method:** The format is taken from `options['format']` (`map`, `newick`, `nexus`, `paf`, `ped`, `phylip`, `twobit` or `zarr`; `fasta` only as the source of `convert`) or detected, first from the path (suffixes `.map`, `.ped`, `.paf`, `.nwk`/`.newick`/`.tre`/`.tree`, `.nex`/`.nexus`/`.nxs`, `.phy`/`.phylip`, `.2bit`, `.fa`/`.fasta`/`.fna`, with a trailing `.gz` ignored; a directory is taken as zarr) and otherwise from the first 64 KiB of content (2bit signature, `#NEXUS`, a leading `(` for Newick, a leading `>` for FASTA, a first line with at least 12 tab-separated columns whose columns 2 to 4 are integers for PAF, an `ntax nchar` header for PHYLIP, four columns ending in an integer for MAP, an even number of at least six columns for PED). The call is then executed by the same engine as the format-specific function (for example `index_nexus`), so the options and the `result` are those of that function. If the format cannot be detected `InputFormatError` (`IO_FORMATS_UNKNOWN_FORMAT`) is raised. A FASTA input with any operation other than `convert` raises `ConfigurationError`.

- **API:** `dnakit.io_formats.index(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a file of any supported format, as accepted by its format-specific function (path, text or bytes). `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Besides the keys below, the options of the format-specific function of the detected format are accepted. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `format`: name of the input format (map, newick, nexus, paf, ped, phylip, twobit, zarr); detected from the suffix/content when omitted

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: exactly the `result` of the format-specific function that was used; its `format` field tells which one (for FASTA input to `convert`, `source_format` is `fasta`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import index

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = index(NEXUS, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "key": "taxon name",
    "record_count": 3,
    "entries": {
      "A": {
        "ordinal": 0,
        "length": 8,
        "residues": 8
      },
      "B": {
        "ordinal": 1,
        "length": 8,
        "residues": 8
      },
      "C d": {
        "ordinal": 2,
        "length": 8,
        "residues": 7
      }
    },
    "trees": {
      "t1": 0
    }
  }
}
```

## 47) index_map

- **Purpose:** Build a lookup index of a PLINK `.map` file: variant IDs with line numbers and byte offsets.

- **Method:** Maps each variant ID to its 1-based line number, byte offset and byte length in the (decompressed) text, together with its chromosome and position, and lists per-chromosome counts, bp ranges and the first line, so a variant can be located without re-parsing. If an ID occurs on several lines, `entries` keeps only its last line (while `record_count` counts every line); such IDs are listed in `duplicate_ids` and a warning is returned.

- **API:** `dnakit.io_formats.index_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.map` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `key` (`"id"`), `record_count`, `entries` (per ID `line`, `offset`, `length`, `chrom`, `pos`) `chromosomes` (per chromosome `variant_count`, `min_pos`, `max_pos`, `first_line`) and `duplicate_ids` (the IDs that occur on more than one line, in order of first occurrence; empty when there are none). When IDs are duplicated, `result` also has `warnings` (a list with one message) and the same message is in the envelope's `warnings`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import index_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
    "X\trs1\t0\t7\n"
)

envelope = index_map(MAP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [
    "1 variant ID(s) occur on several lines (rs1): each entry keeps the last line of its ID, while record_count (5) counts every line."
  ],
  "result": {
    "format": "map",
    "key": "id",
    "record_count": 5,
    "entries": {
      "rs1": {
        "line": 5,
        "offset": 46,
        "length": 10,
        "chrom": "X",
        "pos": 7
      },
      "rs2": {
        "line": 2,
        "offset": 12,
        "length": 13,
        "chrom": "1",
        "pos": 50
      },
      "rs3": {
        "line": 3,
        "offset": 25,
        "length": 11,
        "chrom": "2",
        "pos": 10
      },
      "rs4": {
        "line": 4,
        "offset": 36,
        "length": 10,
        "chrom": "X",
        "pos": 5
      }
    },
    "chromosomes": {
      "1": {
        "variant_count": 2,
        "min_pos": 50,
        "max_pos": 100,
        "first_line": 1
      },
      "2": {
        "variant_count": 1,
        "min_pos": 10,
        "max_pos": 10,
        "first_line": 3
      },
      "X": {
        "variant_count": 2,
        "min_pos": 5,
        "max_pos": 7,
        "first_line": 4
      }
    },
    "duplicate_ids": [
      "rs1"
    ],
    "warnings": [
      "1 variant ID(s) occur on several lines (rs1): each entry keeps the last line of its ID, while record_count (5) counts every line."
    ]
  }
}
```

## 48) index_newick

- **Purpose:** Build a lookup index of a Newick file: trees with byte offsets and the trees containing each tip.

- **Method:** Gives each tree its 0-based index, byte offset and byte length in the (decompressed) text and its tip count, and lists, for every tip name, the indices of the trees that contain it.

- **API:** `dnakit.io_formats.index_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a Newick tree file (one or more `;`-terminated trees), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `key` (`"tree index"`), `record_count`, `trees` (per tree `index`, `offset`, `length`, `tip_count`) and `tips` (tip name to list of tree indices). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import index_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = index_newick(TREES, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "key": "tree index",
    "record_count": 2,
    "trees": [
      {
        "index": 0,
        "offset": 0,
        "length": 30,
        "tip_count": 4
      },
      {
        "index": 1,
        "offset": 31,
        "length": 12,
        "tip_count": 4
      }
    ],
    "tips": {
      "A": [
        0,
        1
      ],
      "B": [
        0,
        1
      ],
      "C": [
        0,
        1
      ],
      "D": [
        0,
        1
      ]
    }
  }
}
```

## 49) index_nexus

- **Purpose:** Build a lookup index of a NEXUS file: taxa and trees.

- **Method:** Maps each taxon name to its row ordinal, sequence length and number of non-missing residues, and each tree name to its position among the trees. There are no byte offsets because NEXUS matrices can be interleaved.

- **API:** `dnakit.io_formats.index_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a NEXUS file (`#NEXUS`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `key` (`"taxon name"`), `record_count`, `entries` (per taxon `ordinal`, `length`, `residues`) and `trees` (tree name to position). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import index_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = index_nexus(NEXUS, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "key": "taxon name",
    "record_count": 3,
    "entries": {
      "A": {
        "ordinal": 0,
        "length": 8,
        "residues": 8
      },
      "B": {
        "ordinal": 1,
        "length": 8,
        "residues": 8
      },
      "C d": {
        "ordinal": 2,
        "length": 8,
        "residues": 7
      }
    },
    "trees": {
      "t1": 0
    }
  }
}
```

## 50) index_paf

- **Purpose:** Build a lookup index of a PAF file: query (or target) names with line numbers, byte offsets and intervals.

- **Method:** With `by='query'` (default) or `by='target'`, maps each sequence name to the list of its alignments, each with the 1-based line number, byte offset and byte length of its line in the (decompressed) text and the aligned interval on that sequence (`start`, `end`, 0-based half-open as in PAF); the list is sorted by start, end and line.

- **API:** `dnakit.io_formats.index_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PAF (Pairwise mApping Format) file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `by`: 'query' (default) or 'target': the name to index alignments by

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `key` (`"query_name"` or `"target_name"`), `record_count` and `entries` (name to list of `{line, offset, length, start, end}`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import index_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = index_paf(PAF, options={"by": "target"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "key": "target_name",
    "record_count": 2,
    "entries": {
      "t1": [
        {
          "line": 1,
          "offset": 0,
          "length": 55,
          "start": 100,
          "end": 190
        }
      ],
      "t2": [
        {
          "line": 2,
          "offset": 55,
          "length": 34,
          "start": 0,
          "end": 50
        }
      ]
    }
  }
}
```

## 51) index_ped

- **Purpose:** Build a lookup index of a PLINK `.ped` file: individuals with line numbers and byte offsets.

- **Method:** Maps `family:individual` to its 1-based line number, byte offset and byte length in the (decompressed) text and lists the individuals of each family.

- **API:** `dnakit.io_formats.index_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.ped` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `max_cells`: int, default 20000000

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `key` (`"family:individual"`), `record_count`, `entries` (per individual `line`, `offset`, `length`) and `families` (family ID to list of individual IDs). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import index_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = index_ped(PED, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "key": "family:individual",
    "record_count": 3,
    "entries": {
      "F1:I1": {
        "line": 1,
        "offset": 0,
        "length": 22
      },
      "F1:I2": {
        "line": 2,
        "offset": 22,
        "length": 22
      },
      "F1:I3": {
        "line": 3,
        "offset": 44,
        "length": 25
      }
    },
    "families": {
      "F1": [
        "I1",
        "I2",
        "I3"
      ]
    }
  }
}
```

## 52) index_phylip

- **Purpose:** Build a lookup index of a PHYLIP file: taxa with row ordinal, length and residue count.

- **Method:** Maps each taxon name to its row ordinal, sequence length and number of non-missing residues. PHYLIP rows are not stored one per line when interleaved, so no byte offsets are given.

- **API:** `dnakit.io_formats.index_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PHYLIP alignment file (header `ntax nchar`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `datatype`: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing characters (gap, ?, N for nucleotides, X for protein)
    - `input_layout`: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
    - `input_strict_names`: bool, default False: input names occupy the first 10 characters of each row (otherwise the first whitespace-delimited word)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `key` (`"taxon name"`), `record_count` and `entries` (per taxon `ordinal`, `length`, `residues`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import index_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = index_phylip(PHYLIP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "key": "taxon name",
    "record_count": 3,
    "entries": {
      "A": {
        "ordinal": 0,
        "length": 8,
        "residues": 8
      },
      "B": {
        "ordinal": 1,
        "length": 8,
        "residues": 8
      },
      "C": {
        "ordinal": 2,
        "length": 8,
        "residues": 7
      }
    }
  }
}
```

## 53) index_twobit

- **Purpose:** Build a lookup index of a 2bit file: sequence names with file offsets.

- **Method:** Reads only the header and the per-sequence record headers: for each sequence the file offset of its record, its length, the number of N blocks and mask blocks and the offset of its packed data.

- **API:** `dnakit.io_formats.index_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a UCSC 2bit genome file, given as a path (`str` or `Path`) or as its raw `bytes` (a `str` that is not an existing path is rejected, because text cannot hold a binary file); both byte orders are accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):


- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `key` (`"sequence name"`), `record_count` and `entries` (per sequence `offset`, `length`, `n_block_count`, `mask_block_count`, `packed_offset`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import base64
import json

from dnakit.io_formats import index_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = index_twobit(TWOBIT, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "key": "sequence name",
    "record_count": 2,
    "entries": {
      "s1": {
        "offset": 30,
        "length": 10,
        "n_block_count": 1,
        "mask_block_count": 1,
        "packed_offset": 62
      },
      "s2": {
        "offset": 65,
        "length": 8,
        "n_block_count": 0,
        "mask_block_count": 0,
        "packed_offset": 81
      }
    }
  }
}
```

## 54) index_zarr

- **Purpose:** Build a lookup index of a Zarr store: array paths with shape and chunk grid.

- **Method:** Lists every array path with its kind, shape, dtype, chunk shape, chunk grid and chunk counts (initialised and total). No chunk data is read.

- **API:** `dnakit.io_formats.index_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the path (`str` or `Path`) of a Zarr directory store (zarr v2 or v3; needs the optional `zarr>=3` and `numpy` packages); text content, bytes and gzip are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):


- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `key` (`"array path"`), `record_count`, `entries` (path to `{shape, dtype, chunks, chunk_grid, nchunks, nchunks_initialized}`), `groups` and `path`. If `zarr` (>= 3) or `numpy` cannot be imported the function returns `status: "not_implemented"` with `required_dependency: "pip install 'zarr>=3' numpy"` and a `message` (`zarr and numpy are required: <import error>`, or `zarr <version> is installed but this adapter needs zarr >= 3.`) instead of a `result`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import index_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = index_zarr("demo.zarr", verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "key": "array path",
    "record_count": 2,
    "entries": {
      "dosage": {
        "shape": [
          4,
          2
        ],
        "dtype": "int64",
        "chunks": [
          4,
          2
        ],
        "chunk_grid": [
          1,
          1
        ],
        "nchunks": 1,
        "nchunks_initialized": 1
      },
      "pos": {
        "shape": [
          4
        ],
        "dtype": "int64",
        "chunks": [
          4
        ],
        "chunk_grid": [
          1
        ],
        "nchunks": 1,
        "nchunks_initialized": 1
      }
    },
    "groups": [],
    "path": "demo.zarr"
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

- **Purpose:** Merge several files of the same supported format into one.

- **Method:** The format is taken from `options['format']` (`map`, `newick`, `nexus`, `paf`, `ped`, `phylip`, `twobit` or `zarr`; `fasta` only as the source of `convert`) or detected, first from the path (suffixes `.map`, `.ped`, `.paf`, `.nwk`/`.newick`/`.tre`/`.tree`, `.nex`/`.nexus`/`.nxs`, `.phy`/`.phylip`, `.2bit`, `.fa`/`.fasta`/`.fna`, with a trailing `.gz` ignored; a directory is taken as zarr) and otherwise from the first 64 KiB of content (2bit signature, `#NEXUS`, a leading `(` for Newick, a leading `>` for FASTA, a first line with at least 12 tab-separated columns whose columns 2 to 4 are integers for PAF, an `ntax nchar` header for PHYLIP, four columns ending in an integer for MAP, an even number of at least six columns for PED). The call is then executed by the same engine as the format-specific function (for example `merge_phylip`), so the options and the `result` are those of that function. If the format cannot be detected `InputFormatError` (`IO_FORMATS_UNKNOWN_FORMAT`) is raised. The input is a non-empty list; the first element decides the format, so all elements must be of that format. A FASTA input with any operation other than `convert` raises `ConfigurationError`.

- **API:** `dnakit.io_formats.merge(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a non-empty list of files of one supported format (or `{'inputs': [...]}`), each as accepted by its format-specific function. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Besides the keys below, the options of the format-specific function of the detected format are accepted. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `format`: name of the input format (map, newick, nexus, paf, ped, phylip, twobit, zarr); detected from the suffix/content when omitted

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: exactly the `result` of the format-specific function that was used; its `format` field tells which one (for FASTA input to `convert`, `source_format` is `fasta`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import merge

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = merge([PHYLIP, PHYLIP], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "record_count": 3,
    "input_count": 2,
    "bytes": 90,
    "sha256": "6ef4966d5f176941e22c4ba195bd108ca859a24e6fcce8f468255ce2bbc51031",
    "text": " 3 16\nA          ACGTACGTACGTACGT\nB          ACGTACGAACGTACGA\nC          ACG-ACGTACG-ACGT\n"
  }
}
```

## 57) merge_map

- **Purpose:** Merge several `.map` files into one PLINK `.map` file.

- **Method:** Variants are concatenated in input order. A variant ID that occurs twice raises `InputFormatError` (`IO_FORMATS_DUPLICATE_ID`) with `duplicates='error'`, or the first occurrence is kept with `duplicates='first'`. With `sort=True` the merged variants are ordered by chromosome (natural order) and bp, then ID.

- **API:** `dnakit.io_formats.merge_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a non-empty list (or `{'inputs': [...]}`) whose items are each a PLINK `.map` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `duplicates`: 'error' (default) or 'first': what to do when a variant ID occurs twice
    - `sort`: bool, default False: sort the merged variants by chromosome and bp

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count` (number of merged inputs), `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import merge_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)
MAP_B = "2\trs3\t0\t10\n3\trs5\t0\t77\n"

envelope = merge_map([MAP, MAP_B], options={"duplicates": "first", "sort": True}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "record_count": 5,
    "input_count": 2,
    "bytes": 57,
    "sha256": "fbdbc7656e124d28d44c918b51a2a367f6266809d601f238fcb2153f6689a483",
    "text": "1\trs2\t0.5\t50\n1\trs1\t0\t100\n2\trs3\t0\t10\n3\trs5\t0\t77\nX\trs4\t0\t5\n"
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

- **Purpose:** Merge several Newick files into one multi-tree Newick file.

- **Method:** The trees of all inputs are concatenated in input order, one normalised Newick line per tree. With `deduplicate=True` a tree whose normalised Newick equals that of an earlier tree is dropped.

- **API:** `dnakit.io_formats.merge_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a non-empty list (or `{'inputs': [...]}`) whose items are each a Newick tree file (one or more `;`-terminated trees), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `deduplicate`: bool, default False: drop trees whose normalised Newick equals an earlier tree

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count` (number of merged inputs), `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import merge_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = merge_newick([TREES, TREES], options={"deduplicate": True}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "record_count": 2,
    "input_count": 2,
    "bytes": 44,
    "sha256": "060f56a376c85eddacd58d52e718831c6b84a379f89394849619c36ed07806e7",
    "text": "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"
  }
}
```

## 60) merge_nexus

- **Purpose:** Merge several NEXUS files into one NEXUS file, concatenating sites per taxon or appending taxa.

- **Method:** `mode='concatenate'` (default) builds a supermatrix: taxa are the union of the inputs in order of first appearance, the sites of same-named taxa are joined, and a taxon absent from an input is filled with the missing character for that input's site count. `mode='append'` adds the taxa of all inputs (duplicate names raise `InputFormatError`; all inputs must have the same number of sites). The datatypes must agree; trees of all inputs are carried over (renamed on collision), without being pruned or extended to the new taxa.

- **API:** `dnakit.io_formats.merge_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a non-empty list (or `{'inputs': [...]}`) whose items are each a NEXUS file (`#NEXUS`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `mode`: 'concatenate' (default: join the sites of the same taxa into a supermatrix, filling taxa absent from an input with the missing character) or 'append' (add taxa; site counts must be equal)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count` (number of merged inputs), `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import merge_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""
NEXUS_B = "#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=2 NCHAR=3;\n  FORMAT DATATYPE=DNA MISSING=? GAP=-;\n  MATRIX\n    A TTT\n    D GGG\n  ;\nEND;\n"

envelope = merge_nexus([NEXUS, NEXUS_B], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "record_count": 5,
    "input_count": 2,
    "bytes": 268,
    "sha256": "18d97702dbe400e22c8bcaf861e10a7269a6320d60863079fe9414e00736a258",
    "text": "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=4 NCHAR=11;\n    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n    MATRIX\n    A      ACGTACGTTTT\n    B      ACGTACGA???\n    'C d'  ACG-ACGT???\n    D      ????????GGG\n    ;\nEND;\n\nBEGIN TREES;\n    TREE t1 = [&R] ((A:1,B:1):1,'C d':2);\nEND;\n"
  }
}
```

## 61) merge_paf

- **Purpose:** Merge several PAF files into one PAF file.

- **Method:** Alignments are concatenated in input order. With `deduplicate=True` a line identical to an earlier one (all fields and tags) is dropped. `sort_by` (a key of `sort_paf`) sorts the merged alignments afterwards.

- **API:** `dnakit.io_formats.merge_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a non-empty list (or `{'inputs': [...]}`) whose items are each a PAF (Pairwise mApping Format) file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `deduplicate`: bool, default False: drop alignment lines identical to an earlier one
    - `sort_by`: optional sort key as for sort_paf, applied to the merged records

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count` (number of merged inputs), `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import merge_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = merge_paf([PAF, PAF], options={"deduplicate": True, "sort_by": "mapq"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "record_count": 2,
    "input_count": 2,
    "bytes": 89,
    "sha256": "5ced75fd261d6ae80944159151cb702eb1c8258e15679a6fba89f91a2282bf99",
    "text": "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\nq1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
  }
}
```

## 62) merge_ped

- **Purpose:** Merge several `.ped` files with the same variants into one PLINK `.ped` file.

- **Method:** The rows are concatenated in input order. All files must have the same number of variants, otherwise `InputFormatError`. A (family, individual) pair that occurs twice raises `InputFormatError` (`IO_FORMATS_DUPLICATE_ID`) with `duplicates='error'`, or the first row is kept with `duplicates='first'`. The function does not check that the variants are the same variants in the same order; that is the caller's responsibility.

- **API:** `dnakit.io_formats.merge_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a non-empty list (or `{'inputs': [...]}`) whose items are each a PLINK `.ped` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `duplicates`: 'error' (default) or 'first': what to do when a (family, individual) pair occurs twice
    - `max_cells`: int, default 20000000

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count` (number of merged inputs), `variant_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import merge_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)
PED_B = "F2 J1 0 0 2 1 C C T T\n"

envelope = merge_ped([PED, PED_B], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "record_count": 4,
    "variant_count": 2,
    "input_count": 2,
    "bytes": 91,
    "sha256": "32ddebc1aec0e9da60a191a2b5a8775a2b4c7151cf4d3be89f05c343daca7a90",
    "text": "F1\tI1\t0\t0\t1\t2\tA\tA\tC\tG\nF1\tI2\t0\t0\t2\t1\tA\tG\t0\t0\nF1\tI3\tI1\tI2\t1\t-9\tG\tG\tC\tC\nF2\tJ1\t0\t0\t2\t1\tC\tC\tT\tT\n"
  }
}
```

## 63) merge_phylip

- **Purpose:** Merge several PHYLIP alignments into one PHYLIP file, concatenating sites per taxon or appending taxa.

- **Method:** `mode='concatenate'` (default) builds a supermatrix: taxa are the union of the inputs in order of first appearance, the sites of same-named taxa are joined, and a taxon absent from an input is filled with `?` for that input's site count. `mode='append'` adds the taxa of all inputs (duplicate names raise `InputFormatError`; all inputs must have the same number of sites). The datatypes must agree.

- **API:** `dnakit.io_formats.merge_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a non-empty list (or `{'inputs': [...]}`) whose items are each a PHYLIP alignment file (header `ntax nchar`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `datatype`: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing characters (gap, ?, N for nucleotides, X for protein)
    - `input_layout`: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
    - `input_strict_names`: bool, default False: input names occupy the first 10 characters of each row (otherwise the first whitespace-delimited word)
    - `mode`: 'concatenate' (default: join the sites of the same taxa into a supermatrix, filling taxa absent from an input with the missing character) or 'append' (add taxa; site counts must be equal)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count` (number of merged inputs), `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import merge_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"
PHYLIP_B = " 2 3\nA  TTT\nD  GGG\n"

envelope = merge_phylip([PHYLIP, PHYLIP_B], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "record_count": 4,
    "input_count": 2,
    "bytes": 98,
    "sha256": "613e33924aea618ca1521ba246371fbe922448952251c676c6d39c4a77ead034",
    "text": " 4 11\nA          ACGTACGTTTT\nB          ACGTACGA???\nC          ACG-ACGT???\nD          ????????GGG\n"
  }
}
```

## 64) merge_twobit

- **Purpose:** Merge several 2bit files into one 2bit file.

- **Method:** Concatenates the sequences of all inputs, in input order, into one version-0 file; the packed records are copied without re-encoding. A name that occurs twice raises `InputFormatError` (`IO_FORMATS_DUPLICATE_ID`) with `duplicates='error'`, or the first sequence is kept with `duplicates='first'`. `byte_order` selects the byte order of the result.

- **API:** `dnakit.io_formats.merge_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a non-empty list (or `{'inputs': [...]}`) whose items are each a UCSC 2bit genome file, given as a path (`str` or `Path`) or as its raw `bytes` (a `str` that is not an existing path is rejected, because text cannot hold a binary file); both byte orders are accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `duplicates`: 'error' (default) or 'first': what to do when a name occurs twice
    - `byte_order`: 'little' (default) or 'big': byte order of the written 2bit file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count` (number of merged inputs), `sequence_count`, `bytes` and `sha256` of the serialised output, and the content as `data_base64`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `data_base64`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import base64
import json

from dnakit.io_formats import merge_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = merge_twobit([TWOBIT, TWOBIT], options={"duplicates": "first"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "record_count": 2,
    "sequence_count": 2,
    "input_count": 2,
    "bytes": 83,
    "sha256": "12ae7a3592d786f9ae950efd4bc62172c1006c9fedd5815f798de219e356f5a1",
    "data_base64": "QydBGgAAAAACAAAAAAAAAAJzMR4AAAACczJBAAAACgAAAAEAAAAEAAAAAgAAAAEAAAAGAAAABAAAAAAAAACcCcAIAAAAAAAAAAAAAAAAAAAA9aA="
  }
}
```

## 65) merge_zarr

- **Purpose:** Merge several Zarr stores by concatenating arrays of the same path into one store.

- **Method:** Every input store must contain the same array paths; the arrays of the same path are concatenated along `axis` (default 0) in input order, shapes must agree except along that axis, and dtypes are promoted with `numpy.result_type`. The result is a new store at `output` (required; `overwrite=True` replaces an existing one). Merging more than `max_elements` elements is refused.

- **API:** `dnakit.io_formats.merge_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a non-empty list (or `{'inputs': [...]}`) whose items are each the path (`str` or `Path`) of a Zarr directory store (zarr v2 or v3; needs the optional `zarr>=3` and `numpy` packages); text content, bytes and gzip are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `output`: directory path of the new zarr store (required)
    - `overwrite`: bool, default False: replace an existing zarr store at output
    - `axis`: int, default 0: axis to concatenate along
    - `max_elements`: int, default 5000000: refuse to materialise more array elements than this

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count` (number of arrays written), `input_count`, `axis` and `output_path` (no content is returned). If `zarr` (>= 3) or `numpy` cannot be imported the function returns `status: "not_implemented"` with `required_dependency: "pip install 'zarr>=3' numpy"` and a `message` (`zarr and numpy are required: <import error>`, or `zarr <version> is installed but this adapter needs zarr >= 3.`) instead of a `result`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import merge_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = merge_zarr(["demo.zarr", "demo.zarr"], options={"output": "merged.zarr", "overwrite": True}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "record_count": 2,
    "input_count": 2,
    "axis": 0,
    "output_path": "merged.zarr"
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

- **Purpose:** Parse a PLINK `.map` file into structured variants.

- **Method:** Each non-blank line must have four whitespace-separated columns: chromosome, variant ID, genetic distance in cM (a finite number) and base-pair position (an integer >= 0). Blank lines are skipped, a malformed line raises `InputFormatError`, and records keep the file order. At most `max_records` variants are returned; use `stream_map` for larger files.

- **API:** `dnakit.io_formats.read_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.map` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `max_records`: int, default 200000: raise instead of returning more records (use stream_* to page)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count` and `records`; each record is `{chrom, id, cm, pos}` (`cm` is a float, `pos` an integer). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import read_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = read_map(MAP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "record_count": 4,
    "records": [
      {
        "chrom": "1",
        "id": "rs1",
        "cm": 0.0,
        "pos": 100
      },
      {
        "chrom": "1",
        "id": "rs2",
        "cm": 0.5,
        "pos": 50
      },
      {
        "chrom": "2",
        "id": "rs3",
        "cm": 0.0,
        "pos": 10
      },
      {
        "chrom": "X",
        "id": "rs4",
        "cm": 0.0,
        "pos": 5
      }
    ]
  }
}
```

## 68) read_newick

- **Purpose:** Parse a Newick tree file into structured trees.

- **Method:** Every `;`-terminated tree is parsed (labels, branch lengths, `[...]` comments) into a node tree; the record carries the re-serialised normalised Newick string, the tip names and a flat node table (`id`, `parent`, `name`, `length`, in preorder, the root has `parent: null`). A syntax error raises `InputFormatError`; trees are limited to `max_nodes` nodes and the file to `max_records` trees.

- **API:** `dnakit.io_formats.read_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a Newick tree file (one or more `;`-terminated trees), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `max_records`: int, default 200000: raise instead of returning more records (use stream_* to page)
    - `max_nodes`: int, default 2000000: maximum nodes per tree

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count` and `records`; each record is `{index, newick, tips, nodes}` with the normalised Newick string, the tip names and the flat node table `{id, parent, name, length}` in preorder. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import read_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = read_newick(TREES, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "record_count": 2,
    "records": [
      {
        "index": 0,
        "newick": "((A:1,B:2):0.5,(C:1,D:1):0.5);",
        "tips": [
          "A",
          "B",
          "C",
          "D"
        ],
        "nodes": [
          {
            "id": 0,
            "parent": null,
            "name": null,
            "length": null
          },
          {
            "id": 1,
            "parent": 0,
            "name": null,
            "length": 0.5
          },
          {
            "id": 2,
            "parent": 1,
            "name": "A",
            "length": 1.0
          },
          {
            "id": 3,
            "parent": 1,
            "name": "B",
            "length": 2.0
          },
          {
            "id": 4,
            "parent": 0,
            "name": null,
            "length": 0.5
          },
          {
            "id": 5,
            "parent": 4,
            "name": "C",
            "length": 1.0
          },
          {
            "id": 6,
            "parent": 4,
            "name": "D",
            "length": 1.0
          }
        ]
      },
      {
        "index": 1,
        "newick": "(A,B,(C,D));",
        "tips": [
          "A",
          "B",
          "C",
          "D"
        ],
        "nodes": [
          {
            "id": 0,
            "parent": null,
            "name": null,
            "length": null
          },
          {
            "id": 1,
            "parent": 0,
            "name": "A",
            "length": null
          },
          {
            "id": 2,
            "parent": 0,
            "name": "B",
            "length": null
          },
          {
            "id": 3,
            "parent": 0,
            "name": null,
            "length": null
          },
          {
            "id": 4,
            "parent": 3,
            "name": "C",
            "length": null
          },
          {
            "id": 5,
            "parent": 3,
            "name": "D",
            "length": null
          }
        ]
      }
    ]
  }
}
```

## 69) read_nexus

- **Purpose:** Parse a NEXUS file into structured sequences and trees.

- **Method:** Parses the TAXA, DATA/CHARACTERS and TREES blocks (other blocks are ignored): taxon labels, the matrix (sequential or interleaved, `matchchar` expanded), the declared datatype, missing and gap symbols, and the trees (TRANSLATE tables applied, `[&R]`/`[&U]` kept as the `rooted` flag). Comments in `[]` are dropped. Quoted taxon names (`'C d'`) are unquoted.

- **API:** `dnakit.io_formats.read_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a NEXUS file (`#NEXUS`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `max_records`: int, default 200000: raise instead of returning more records (use stream_* to page)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `taxon_count`, `site_count`, `datatype`, `taxa`, `tree_count`, `missing`, `gap`, `layout` and `records` (sequence records `{kind: "sequence", name, sequence}` followed by tree records `{kind: "tree", name, rooted, newick}`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import read_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = read_nexus(NEXUS, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "record_count": 4,
    "taxon_count": 3,
    "site_count": 8,
    "datatype": "dna",
    "taxa": [
      "A",
      "B",
      "C d"
    ],
    "tree_count": 1,
    "missing": "?",
    "gap": "-",
    "layout": "sequential",
    "records": [
      {
        "kind": "sequence",
        "name": "A",
        "sequence": "ACGTACGT"
      },
      {
        "kind": "sequence",
        "name": "B",
        "sequence": "ACGTACGA"
      },
      {
        "kind": "sequence",
        "name": "C d",
        "sequence": "ACG-ACGT"
      },
      {
        "kind": "tree",
        "name": "t1",
        "rooted": true,
        "newick": "((A:1,B:1):1,'C d':2);"
      }
    ]
  }
}
```

## 70) read_paf

- **Purpose:** Parse a PAF alignment file into structured alignments.

- **Method:** Each non-blank line is split on tabs into the 12 mandatory columns (query name/length/start/end, strand, target name/length/start/end, matching bases, block length, mapping quality); further columns are optional SAM-like tags `TAG:TYPE:VALUE`, returned as `tags` with the declared type kept next to the value. A malformed line raises `InputFormatError`. At most `max_records` alignments are returned.

- **API:** `dnakit.io_formats.read_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PAF (Pairwise mApping Format) file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `max_records`: int, default 200000: raise instead of returning more records (use stream_* to page)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count` and `records`; each record has the 12 PAF fields plus `tags` (tag to `{type, value}`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import read_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = read_paf(PAF, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "record_count": 2,
    "records": [
      {
        "query_name": "q1",
        "query_length": 100,
        "query_start": 0,
        "query_end": 90,
        "strand": "+",
        "target_name": "t1",
        "target_length": 1000,
        "target_start": 100,
        "target_end": 190,
        "matches": 80,
        "block_length": 90,
        "mapq": 60,
        "tags": {
          "tp": {
            "type": "A",
            "value": "P"
          },
          "cg": {
            "type": "Z",
            "value": "90M"
          }
        }
      },
      {
        "query_name": "q2",
        "query_length": 50,
        "query_start": 0,
        "query_end": 50,
        "strand": "-",
        "target_name": "t2",
        "target_length": 500,
        "target_start": 0,
        "target_end": 50,
        "matches": 40,
        "block_length": 50,
        "mapq": 10,
        "tags": {}
      }
    ]
  }
}
```

## 71) read_ped

- **Purpose:** Parse a PLINK `.ped` file into structured individuals.

- **Method:** Each row has `6 + 2k` whitespace-separated columns (family, individual, father, mother, sex, phenotype, then two alleles per variant) with the same `k` on every row. `sex` is read as an integer; `pheno` is kept as the original string. The genotype matrix is bounded by `max_cells` (individuals x variants) and the number of rows by `max_records`.

- **API:** `dnakit.io_formats.read_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.ped` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `max_records`: int, default 200000: raise instead of returning more records (use stream_* to page)
    - `max_cells`: int, default 20000000: refuse files with more genotype cells (individuals x variants)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `variant_count` and `records`; each record is `{fid, iid, pat, mat, sex, pheno, genotypes}` with `genotypes` a list of `[allele1, allele2]` pairs (`sex` is an integer, `pheno` a string; allele `0` marks a missing call). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import read_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = read_ped(PED, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "record_count": 3,
    "variant_count": 2,
    "records": [
      {
        "fid": "F1",
        "iid": "I1",
        "pat": "0",
        "mat": "0",
        "sex": 1,
        "pheno": "2",
        "genotypes": [
          [
            "A",
            "A"
          ],
          [
            "C",
            "G"
          ]
        ]
      },
      {
        "fid": "F1",
        "iid": "I2",
        "pat": "0",
        "mat": "0",
        "sex": 2,
        "pheno": "1",
        "genotypes": [
          [
            "A",
            "G"
          ],
          [
            "0",
            "0"
          ]
        ]
      },
      {
        "fid": "F1",
        "iid": "I3",
        "pat": "I1",
        "mat": "I2",
        "sex": 1,
        "pheno": "-9",
        "genotypes": [
          [
            "G",
            "G"
          ],
          [
            "C",
            "C"
          ]
        ]
      }
    ]
  }
}
```

## 72) read_phylip

- **Purpose:** Parse a PHYLIP alignment file into structured sequences.

- **Method:** Reads relaxed PHYLIP (the name is the first whitespace-delimited word; residue blocks may contain spaces) or strict PHYLIP (`input_strict_names=True`: the name is the first 10 characters); sequential and interleaved layouts are detected automatically unless `input_layout` is given. The header `ntax nchar` must match the data. `datatype` selects the alphabet used for validation.

- **API:** `dnakit.io_formats.read_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PHYLIP alignment file (header `ntax nchar`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `max_records`: int, default 200000: raise instead of returning more records (use stream_* to page)
    - `datatype`: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing characters (gap, ?, N for nucleotides, X for protein)
    - `input_layout`: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
    - `input_strict_names`: bool, default False: input names occupy the first 10 characters of each row (otherwise the first whitespace-delimited word)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `taxon_count`, `site_count`, `datatype` and `records` (each record is `{name, sequence}`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import read_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = read_phylip(PHYLIP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "record_count": 3,
    "taxon_count": 3,
    "site_count": 8,
    "datatype": "dna",
    "records": [
      {
        "name": "A",
        "sequence": "ACGTACGT"
      },
      {
        "name": "B",
        "sequence": "ACGTACGA"
      },
      {
        "name": "C",
        "sequence": "ACG-ACGT"
      }
    ]
  }
}
```

## 73) read_twobit

- **Purpose:** Parse a UCSC 2bit file into sequences, or into a region of one sequence.

- **Method:** Decodes the packed 2-bit bases to text: N blocks become `N` and mask blocks become lower case (`soft_mask=False` gives upper case). `names` selects sequences (default all); `start`/`end` give the 0-based half-open region `[start, end)` of exactly one named sequence, otherwise `ConfigurationError`. Decoding more than `max_bases` bases is refused.

- **API:** `dnakit.io_formats.read_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a UCSC 2bit genome file, given as a path (`str` or `Path`) or as its raw `bytes` (a `str` that is not an existing path is rejected, because text cannot hold a binary file); both byte orders are accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `names`: list of sequence names to read (default all)
    - `start`: int: 0-based start (single name only)
    - `end`: int: 0-based exclusive end (single name only)
    - `soft_mask`: bool, default True: lower-case masked bases
    - `max_bases`: int, default 50000000: refuse to decode more bases

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count` and `records`; each record is `{name, start, end, length, sequence}` where `length` is the full sequence length and `[start, end)` the decoded region. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import base64
import json

from dnakit.io_formats import read_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = read_twobit(TWOBIT, options={"names": ["s1"], "start": 2, "end": 8}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "records": [
      {
        "name": "s1",
        "start": 2,
        "end": 8,
        "length": 10,
        "sequence": "GTNNac"
      }
    ],
    "record_count": 1
  }
}
```

## 74) read_zarr

- **Purpose:** Read a Zarr directory store into structured arrays.

- **Method:** Opens the store read-only and returns, per array, shape, dtype, chunk shape, chunk grid, chunk counts, compressors, zarr format, uncompressed size and attributes and, unless `include_data=False`, the values as nested lists. `array` selects one array (default all); `start`/`stop` select rows along axis 0 of a single array. Reading more than `max_elements` elements is refused.

- **API:** `dnakit.io_formats.read_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the path (`str` or `Path`) of a Zarr directory store (zarr v2 or v3; needs the optional `zarr>=3` and `numpy` packages); text content, bytes and gzip are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `array`: array path to read (default: all arrays)
    - `include_data`: bool, default True: return element values
    - `start`: int: first row along axis 0 (single array)
    - `stop`: int: end row (exclusive)
    - `max_elements`: int, default 5000000: refuse to materialise more array elements than this

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `path`, `arrays` (path to `{shape, dtype, chunks, chunk_grid, nchunks, nchunks_initialized, compressors, zarr_format, nbytes, attrs, data}`; `data` is absent with `include_data=False`), `groups`, `attrs`, `record_count` and `records` (the same metadata without `data`, with `path`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import read_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = read_zarr("demo.zarr", options={"array": "dosage", "start": 1, "stop": 3}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "path": "demo.zarr",
    "arrays": {
      "dosage": {
        "shape": [
          4,
          2
        ],
        "dtype": "int64",
        "chunks": [
          4,
          2
        ],
        "chunk_grid": [
          1,
          1
        ],
        "nchunks": 1,
        "nchunks_initialized": 1,
        "compressors": [
          "ZstdCodec"
        ],
        "zarr_format": 3,
        "nbytes": 64,
        "attrs": {},
        "data": [
          [
            1,
            1
          ],
          [
            0,
            0
          ]
        ]
      }
    },
    "groups": [],
    "attrs": {},
    "record_count": 1,
    "records": [
      {
        "path": "dosage",
        "shape": [
          4,
          2
        ],
        "dtype": "int64",
        "chunks": [
          4,
          2
        ],
        "chunk_grid": [
          1,
          1
        ],
        "nchunks": 1,
        "nchunks_initialized": 1,
        "compressors": [
          "ZstdCodec"
        ],
        "zarr_format": 3,
        "nbytes": 64,
        "attrs": {}
      }
    ]
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

- **Purpose:** Download a file over HTTP(S), resuming an interrupted transfer with a Range request.

- **Method:** The URL is fetched into `<output>.part`. If that part file already holds bytes, they are kept and the rest is requested with `Range: bytes=<n>-`, guarded by `If-Range` with the ETag (or Last-Modified) saved in `<output>.part.meta`. A server answering `200` (it ignored Range, or the resource changed) restarts the download from byte 0 (`restarted: true`); `416` with a matching total size means the part file is already complete. A `206` reply must start at the requested offset. When the received size equals the announced size and, if `expected_sha256` is given, the digest matches (a mismatch deletes the part file and raises `DownloadError`), the part file is renamed to `output`. A transfer cut short keeps the `.part` file and raises `DownloadError` (`DOWNLOAD_LENGTH_MISMATCH` or `DOWNLOAD_NETWORK_ERROR`), so calling the function again resumes it. Only `http` and `https` URLs are accepted (also across redirects); more than `max_bytes` raises `DownloadError`. This is the only function of the page that uses the network.

- **API:** `dnakit.io_formats.resume_download(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the URL (`http://` or `https://`) as a string, or a mapping `{'url': ..., 'output': ...}`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `url`: http(s) URL (alternative to passing it as input)
    - `output`: destination file path (required unless input is a mapping with 'output')
    - `overwrite`: bool, default False: replace an existing destination file
    - `chunk_size`: int, default 65536: bytes read per step
    - `timeout`: number of seconds, default 30: socket timeout
    - `max_bytes`: int, default 1073741824: refuse to download more than this
    - `expected_sha256`: hex digest the finished file must have (a mismatch deletes the partial file)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `output_path`, `bytes` (size of the finished file), `sha256`, `restarted` (the download had to start again from byte 0), `complete` (always `true` on success) and `resumed_from` (byte offset the transfer continued from, 0 for a fresh download). HTTP or network failures raise `DownloadError`; invalid options raise `ConfigurationError`; a URL that is not http(s) raises `InputFormatError`.

- **Example code:**

```python
import hashlib
import http.server
import json
import threading
from pathlib import Path

from dnakit.io_formats import resume_download

PAYLOAD = bytes(range(256)) * 40  # the remote file: 10240 bytes


class RangeHandler(http.server.BaseHTTPRequestHandler):
    """A tiny local server that honours 'Range: bytes=<n>-' requests."""

    def log_message(self, *args):
        pass

    def do_GET(self):
        start = 0
        header = self.headers.get("Range")
        if header:
            start = int(header.split("=")[1].rstrip("-"))
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(PAYLOAD) - 1}/{len(PAYLOAD)}")
        else:
            self.send_response(200)
        self.send_header("Content-Length", str(len(PAYLOAD) - start))
        self.send_header("ETag", '"v1"')
        self.end_headers()
        self.wfile.write(PAYLOAD[start:])


server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), RangeHandler)
threading.Thread(target=server.serve_forever, daemon=True).start()
url = f"http://127.0.0.1:{server.server_address[1]}/file.bin"

# an interrupted earlier transfer: 3000 bytes are already stored in out.bin.part
Path("out.bin.part").write_bytes(PAYLOAD[:3000])
Path("out.bin.part.meta").write_text(json.dumps({"etag": '"v1"'}))

envelope = resume_download(
    url,
    options={"output": "out.bin", "expected_sha256": hashlib.sha256(PAYLOAD).hexdigest()},
    verbose=False,
)
server.shutdown()
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
assert Path("out.bin").read_bytes() == PAYLOAD  # the finished file equals the remote file
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "output_path": "out.bin",
    "bytes": 10240,
    "sha256": "e96760a87768717bcebcfd25ddc7d46b4dbc95a4b0014def080c08539f7d90d0",
    "restarted": false,
    "complete": true,
    "resumed_from": 3000
  }
}
```

## 77) sort

- **Purpose:** Sort a file of any supported format with the ordering of its format-specific sort function.

- **Method:** The format is taken from `options['format']` (`map`, `newick`, `nexus`, `paf`, `ped`, `phylip`, `twobit` or `zarr`; `fasta` only as the source of `convert`) or detected, first from the path (suffixes `.map`, `.ped`, `.paf`, `.nwk`/`.newick`/`.tre`/`.tree`, `.nex`/`.nexus`/`.nxs`, `.phy`/`.phylip`, `.2bit`, `.fa`/`.fasta`/`.fna`, with a trailing `.gz` ignored; a directory is taken as zarr) and otherwise from the first 64 KiB of content (2bit signature, `#NEXUS`, a leading `(` for Newick, a leading `>` for FASTA, a first line with at least 12 tab-separated columns whose columns 2 to 4 are integers for PAF, an `ntax nchar` header for PHYLIP, four columns ending in an integer for MAP, an even number of at least six columns for PED). The call is then executed by the same engine as the format-specific function (for example `sort_map`), so the options and the `result` are those of that function. If the format cannot be detected `InputFormatError` (`IO_FORMATS_UNKNOWN_FORMAT`) is raised. A FASTA input with any operation other than `convert` raises `ConfigurationError`.

- **API:** `dnakit.io_formats.sort(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a file of any supported format, as accepted by its format-specific function (path, text or bytes). `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Besides the keys below, the options of the format-specific function of the detected format are accepted. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `format`: name of the input format (map, newick, nexus, paf, ped, phylip, twobit, zarr); detected from the suffix/content when omitted

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: exactly the `result` of the format-specific function that was used; its `format` field tells which one (for FASTA input to `convert`, `source_format` is `fasta`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import sort

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = sort(MAP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "record_count": 4,
    "input_count": 4,
    "output_count": 4,
    "bytes": 46,
    "sha256": "4299c54d1a304472628eaa9018e807d127e591bc315c46a7272a56910661277b",
    "text": "1\trs2\t0.5\t50\n1\trs1\t0\t100\n2\trs3\t0\t10\nX\trs4\t0\t5\n"
  }
}
```

## 78) sort_map

- **Purpose:** Sort the variants of a PLINK `.map` file by chromosome and position (or by ID or genetic distance).

- **Method:** `key='position'` (default) orders by chromosome in natural order (1..22, X, Y, XY, M/MT, then other names alphabetically), then bp, then ID; `key='id'` orders by ID, then chromosome and bp; `key='cm'` orders by genetic distance, then chromosome and bp. `reverse=True` reverses the order.

- **API:** `dnakit.io_formats.sort_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.map` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `key`: 'position' (default: chromosome in natural order, then bp, then ID), 'id' or 'cm'
    - `reverse`: bool, default False

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import sort_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = sort_map(MAP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "record_count": 4,
    "input_count": 4,
    "output_count": 4,
    "bytes": 46,
    "sha256": "4299c54d1a304472628eaa9018e807d127e591bc315c46a7272a56910661277b",
    "text": "1\trs2\t0.5\t50\n1\trs1\t0\t100\n2\trs3\t0\t10\nX\trs4\t0\t5\n"
  }
}
```

## 79) sort_newick

- **Purpose:** Sort the clades of each tree (ladderize or alphabetical order) of a Newick file.

- **Method:** The children of every node are reordered, the tree itself is unchanged: `order='ascending'` (default) puts clades with fewer tips first (ladderize), `'descending'` the reverse, `'alphabetical'` orders clades by their smallest tip name; ties are broken by that name. The order of the trees in the file is kept.

- **API:** `dnakit.io_formats.sort_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a Newick tree file (one or more `;`-terminated trees), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `order`: 'ascending' (default: clades with fewer tips first), 'descending' or 'alphabetical' (by the smallest tip name of each clade); ties are broken by that name

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import sort_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = sort_newick(TREES, options={"order": "descending"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "record_count": 2,
    "input_count": 2,
    "output_count": 2,
    "bytes": 44,
    "sha256": "654ac0beaf1d726fa4cc69552def8d54c2df37576b814570a46d7560004d7032",
    "text": "((A:1,B:2):0.5,(C:1,D:1):0.5);\n((C,D),A,B);\n"
  }
}
```

## 80) sort_nexus

- **Purpose:** Sort the taxa (matrix rows) of a NEXUS file by name or missing-data fraction.

- **Method:** `key='name'` (default) orders rows by taxon name; `key='missing'` by the fraction of gap/missing/ambiguous characters, then name. `reverse=True` reverses the order. Sites and trees are unchanged.

- **API:** `dnakit.io_formats.sort_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a NEXUS file (`#NEXUS`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `key`: 'name' (default) or 'missing'
    - `reverse`: bool, default False

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import sort_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = sort_nexus(NEXUS, options={"key": "missing", "reverse": True}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "record_count": 4,
    "input_count": 4,
    "output_count": 4,
    "bytes": 235,
    "sha256": "d105b30c92592d6e0bc698e1fd80bfebaadaf5d21acb7a30b952b2b45dd18f83",
    "text": "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=3 NCHAR=8;\n    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n    MATRIX\n    'C d'  ACG-ACGT\n    B      ACGTACGA\n    A      ACGTACGT\n    ;\nEND;\n\nBEGIN TREES;\n    TREE t1 = [&R] ((A:1,B:1):1,'C d':2);\nEND;\n"
  }
}
```

## 81) sort_paf

- **Purpose:** Sort the alignments of a PAF file by target, query, mapping quality, identity, block length or matches.

- **Method:** `key='target'` (default) orders by target name, target start and end; `query` by query name, query start and end; `mapq`, `identity` (matches / block_length), `block_length` and `matches` by that value. `reverse=True` reverses the order.

- **API:** `dnakit.io_formats.sort_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PAF (Pairwise mApping Format) file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `key`: 'target' (default: target name, start, end), 'query' (query name, start, end), 'mapq', 'identity', 'block_length' or 'matches'
    - `reverse`: bool, default False

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import sort_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = sort_paf(PAF, options={"key": "mapq"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "record_count": 2,
    "input_count": 2,
    "output_count": 2,
    "bytes": 89,
    "sha256": "5ced75fd261d6ae80944159151cb702eb1c8258e15679a6fba89f91a2282bf99",
    "text": "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\nq1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
  }
}
```

## 82) sort_ped

- **Purpose:** Sort the individuals of a PLINK `.ped` file by family, individual ID, sex or phenotype.

- **Method:** `key='family'` (default) orders by family ID then individual ID; `individual` by individual ID then family ID; `sex` by the sex code, then family and individual; `phenotype` by the phenotype value as a number, then family and individual (the parser requires a numeric phenotype in every row). IDs are compared as strings. `reverse=True` reverses the order.

- **API:** `dnakit.io_formats.sort_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.ped` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `key`: 'family' (default: family ID then individual ID), 'individual', 'sex' or 'phenotype'
    - `reverse`: bool, default False
    - `max_cells`: int, default 20000000

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `variant_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import sort_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = sort_ped(PED, options={"key": "phenotype"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "record_count": 3,
    "variant_count": 2,
    "input_count": 3,
    "output_count": 3,
    "bytes": 69,
    "sha256": "e96809700ea4ad197aae18b2174944b94d7ec502c4febe6d1549d56d1d78f358",
    "text": "F1\tI3\tI1\tI2\t1\t-9\tG\tG\tC\tC\nF1\tI2\t0\t0\t2\t1\tA\tG\t0\t0\nF1\tI1\t0\t0\t1\t2\tA\tA\tC\tG\n"
  }
}
```

## 83) sort_phylip

- **Purpose:** Sort the taxa (rows) of a PHYLIP file by name or missing-data fraction.

- **Method:** `key='name'` (default) orders rows by taxon name; `key='missing'` by the fraction of gap/missing/ambiguous characters, then name. `reverse=True` reverses the order. The sites are unchanged.

- **API:** `dnakit.io_formats.sort_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PHYLIP alignment file (header `ntax nchar`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `datatype`: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing characters (gap, ?, N for nucleotides, X for protein)
    - `input_layout`: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
    - `input_strict_names`: bool, default False: input names occupy the first 10 characters of each row (otherwise the first whitespace-delimited word)
    - `key`: 'name' (default) or 'missing' (missing fraction, then name)
    - `reverse`: bool, default False

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import sort_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = sort_phylip(PHYLIP, options={"key": "missing", "reverse": True}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "record_count": 3,
    "input_count": 3,
    "output_count": 3,
    "bytes": 65,
    "sha256": "05bc6f3276e92b5bef45beaa26c3ffeb9d2af1c36f3e0ec87386a1152c9f7677",
    "text": " 3 8\nC          ACG-ACGT\nB          ACGTACGA\nA          ACGTACGT\n"
  }
}
```

## 84) sort_twobit

- **Purpose:** Sort the sequences of a 2bit file by name, natural chromosome order or length.

- **Method:** The sequences are reordered in the output index and the packed records are copied bit-exactly. `key='name'` (default) is plain name order; `'natural'` is chromosome order (`chr1`..`chr22`, `X`, `Y`, `M`, then other names alphabetically); `'length'` is by sequence length. `reverse=True` reverses the order.

- **API:** `dnakit.io_formats.sort_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a UCSC 2bit genome file, given as a path (`str` or `Path`) or as its raw `bytes` (a `str` that is not an existing path is rejected, because text cannot hold a binary file); both byte orders are accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output`: path to write (optional)
    - `overwrite`: bool, default False
    - `key`: 'name' (default), 'natural' (chr1..22, X, Y, M order) or 'length'
    - `reverse`: bool, default False
    - `byte_order`: 'little' (default) or 'big': byte order of the written 2bit file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `input_count`, `output_count`, `sequence_count`, `bytes` and `sha256` of the serialised output, and the content as `data_base64`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `data_base64`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import base64
import json

from dnakit.io_formats import sort_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = sort_twobit(TWOBIT, options={"key": "length"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "record_count": 2,
    "sequence_count": 2,
    "input_count": 2,
    "output_count": 2,
    "bytes": 83,
    "sha256": "651601b14986cb13fa6b5342f96e6fbb9ed4474860233892def1156920d74518",
    "data_base64": "QydBGgAAAAACAAAAAAAAAAJzMh4AAAACczEwAAAACAAAAAAAAAAAAAAAAAAAAPWgCgAAAAEAAAAEAAAAAgAAAAEAAAAGAAAABAAAAAAAAACcCcA="
  }
}
```

## 85) sort_zarr

- **Purpose:** Sort the rows (axis 0) of all arrays of a Zarr store by the values of a key array, writing a new store.

- **Method:** The 1-D array named by `by` (required) gives a stable ordering of its values (`reverse=True` reverses the ascending order); axis 0 of every array whose first dimension equals the length of the key array is permuted by that order, and other arrays are copied unchanged. The result is a new store at `output` (required). Sorting more than `max_elements` elements is refused.

- **API:** `dnakit.io_formats.sort_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the path (`str` or `Path`) of a Zarr directory store (zarr v2 or v3; needs the optional `zarr>=3` and `numpy` packages); text content, bytes and gzip are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `output`: directory path of the new zarr store (required)
    - `overwrite`: bool, default False: replace an existing zarr store at output
    - `by`: path of the 1-D key array (required); rows are ordered by its values (stable)
    - `reverse`: bool, default False: reverse the ascending order
    - `max_elements`: int, default 5000000: refuse to materialise more array elements than this

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count` (number of arrays written), `key`, `permuted_arrays`, `copied_unchanged` and `output_path` (no content is returned). If `zarr` (>= 3) or `numpy` cannot be imported the function returns `status: "not_implemented"` with `required_dependency: "pip install 'zarr>=3' numpy"` and a `message` (`zarr and numpy are required: <import error>`, or `zarr <version> is installed but this adapter needs zarr >= 3.`) instead of a `result`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import sort_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = sort_zarr("demo.zarr", options={"output": "sorted.zarr", "overwrite": True, "by": "pos"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "record_count": 2,
    "key": "pos",
    "permuted_arrays": [
      "dosage",
      "pos"
    ],
    "copied_unchanged": [],
    "output_path": "sorted.zarr"
  }
}
```

## 86) split_map

- **Purpose:** Split a PLINK `.map` file per chromosome, or into blocks of N variants.

- **Method:** With `by='chromosome'` (default) one part per chromosome is produced, ordered naturally (1..22, X, Y, XY, M/MT, others) and named by the chromosome. With `by='size'` consecutive blocks of `size` variants are produced, named `0001`, `0002`, ...; `size` is then required. Each part is a complete `.map` file.

- **API:** `dnakit.io_formats.split_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.map` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output_dir`: directory to write the parts to as part_0001_<label><ext> (optional)
    - `overwrite`: bool, default False
    - `by`: 'chromosome' (default; one part per chromosome) or 'size'
    - `size`: int: variants per part when by='size'

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `part_count` and `parts`; every part has `name`, `record_count`, `bytes`, `sha256` and its content (`text`, or `data_base64` for binary formats); with `options['output_dir']` each part is written to `part_0001_<name><ext>`, `part_0002_...` and `output_path` is returned instead of the content. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import split_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = split_map(MAP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "part_count": 3,
    "parts": [
      {
        "name": "1",
        "record_count": 2,
        "bytes": 25,
        "sha256": "40be7a615faec1b5f39695fe44a7f975a678b81bb09a0ac0aa2608389cbefd6e",
        "text": "1\trs1\t0\t100\n1\trs2\t0.5\t50\n"
      },
      {
        "name": "2",
        "record_count": 1,
        "bytes": 11,
        "sha256": "fc5af41f3a73f2c542af0c2b9e7be799905f0ae6459c62f6162daa2088473cbd",
        "text": "2\trs3\t0\t10\n"
      },
      {
        "name": "X",
        "record_count": 1,
        "bytes": 10,
        "sha256": "ac4cc4dc84e656a5e8cf48bab597529c20e7a0af3d275866c00dccc7470660bb",
        "text": "X\trs4\t0\t5\n"
      }
    ]
  }
}
```

## 87) split_newick

- **Purpose:** Split a Newick file into files of N trees.

- **Method:** Consecutive groups of `trees_per_part` trees (default 1) become the parts, named `0001`, `0002`, ...; each part is a complete Newick file.

- **API:** `dnakit.io_formats.split_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a Newick tree file (one or more `;`-terminated trees), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output_dir`: directory to write the parts to as part_0001_<label><ext> (optional)
    - `overwrite`: bool, default False
    - `trees_per_part`: int, default 1: trees in each output part

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `part_count` and `parts`; every part has `name`, `record_count`, `bytes`, `sha256` and its content (`text`, or `data_base64` for binary formats); with `options['output_dir']` each part is written to `part_0001_<name><ext>`, `part_0002_...` and `output_path` is returned instead of the content. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import split_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = split_newick(TREES, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "part_count": 2,
    "parts": [
      {
        "name": "0001",
        "record_count": 1,
        "bytes": 31,
        "sha256": "2443c552d427d6bcde12f7e9de1eea8a374ca83fae68f6b5be7f2cdb9de0b008",
        "text": "((A:1,B:2):0.5,(C:1,D:1):0.5);\n"
      },
      {
        "name": "0002",
        "record_count": 1,
        "bytes": 13,
        "sha256": "2cead55aa9365339f55531d03a9973092b4c1d2a66503e615296d53c826612b6",
        "text": "(A,B,(C,D));\n"
      }
    ]
  }
}
```

## 88) split_nexus

- **Purpose:** Split a NEXUS file by site blocks, taxon groups, named site partitions or tree groups.

- **Method:** Give exactly one of `site_size` (consecutive blocks of sites, parts named `sites_<first>-<last>`), `taxa_size` (consecutive groups of taxa, `taxa_0001`, ...), `partitions` (label to `start-end[,start-end...]`, 1-based inclusive site ranges, one part per label) or `trees_per_part` (groups of trees, `trees_0001`, ...; the file must contain trees). Alignment parts contain the matrix only; tree parts contain a TAXA block and the trees. Each part is a complete NEXUS file.

- **API:** `dnakit.io_formats.split_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a NEXUS file (`#NEXUS`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output_dir`: directory to write the parts to as part_0001_<label><ext> (optional)
    - `overwrite`: bool, default False
    - `site_size`: int: consecutive blocks of this many sites
    - `taxa_size`: int: consecutive groups of this many taxa
    - `partitions`: mapping label -> 'start-end[,start-end...]' (1-based, inclusive site ranges)
    - `trees_per_part`: int: trees per part (documents without a matrix); give exactly one of the four split options

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `part_count` and `parts`; every part has `name`, `record_count`, `bytes`, `sha256` and its content (`text`, or `data_base64` for binary formats); with `options['output_dir']` each part is written to `part_0001_<name><ext>`, `part_0002_...` and `output_path` is returned instead of the content. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import split_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = split_nexus(NEXUS, options={"site_size": 4}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "part_count": 2,
    "parts": [
      {
        "name": "sites_1-4",
        "record_count": 3,
        "bytes": 162,
        "sha256": "f730e6e646487331208af5613fdcf90ac1bcd885e831dc963b4c860a4132b190",
        "text": "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=3 NCHAR=4;\n    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n    MATRIX\n    A      ACGT\n    B      ACGT\n    'C d'  ACG-\n    ;\nEND;\n"
      },
      {
        "name": "sites_5-8",
        "record_count": 3,
        "bytes": 162,
        "sha256": "4f5f5fcb316930bf01371cd7d93e8f613a41f263e25d7d12b1a2522983797e59",
        "text": "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=3 NCHAR=4;\n    FORMAT DATATYPE=DNA MISSING=? GAP=-;\n    MATRIX\n    A      ACGT\n    B      ACGA\n    'C d'  ACGT\n    ;\nEND;\n"
      }
    ]
  }
}
```

## 89) split_paf

- **Purpose:** Split a PAF file per query, target or strand, or into blocks of N alignments.

- **Method:** `by='query'` (default), `'target'` or `'strand'` makes one part per distinct value, in order of first appearance and named by it; `by='size'` makes consecutive blocks of `size` alignments, named `0001`, `0002`, ... (`size` is then required). Each part is a complete PAF file.

- **API:** `dnakit.io_formats.split_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PAF (Pairwise mApping Format) file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output_dir`: directory to write the parts to as part_0001_<label><ext> (optional)
    - `overwrite`: bool, default False
    - `by`: 'query' (default), 'target', 'strand' or 'size'
    - `size`: int: alignments per part when by='size'

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `part_count` and `parts`; every part has `name`, `record_count`, `bytes`, `sha256` and its content (`text`, or `data_base64` for binary formats); with `options['output_dir']` each part is written to `part_0001_<name><ext>`, `part_0002_...` and `output_path` is returned instead of the content. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import split_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = split_paf(PAF, options={"by": "strand"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "part_count": 2,
    "parts": [
      {
        "name": "+",
        "record_count": 1,
        "bytes": 55,
        "sha256": "5f960e3807baf139ee66a38539d822ab8784cde4fe8a933cdbb7c1d941d2231b",
        "text": "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
      },
      {
        "name": "-",
        "record_count": 1,
        "bytes": 34,
        "sha256": "5bec66460b177fb4f757cf63b7dc30daee9f94c45e100137e47512b6dc9b5bd3",
        "text": "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
      }
    ]
  }
}
```

## 90) split_ped

- **Purpose:** Split a PLINK `.ped` file per family, or into blocks of N individuals.

- **Method:** With `by='family'` (default) one part per family is produced, in order of first appearance and named by the family ID. With `by='size'` consecutive blocks of `size` individuals are produced, named `0001`, `0002`, ...; `size` is then required. Every part keeps all variant columns and is a complete `.ped` file.

- **API:** `dnakit.io_formats.split_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.ped` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output_dir`: directory to write the parts to as part_0001_<label><ext> (optional)
    - `overwrite`: bool, default False
    - `by`: 'family' (default; one part per family) or 'size'
    - `size`: int: individuals per part when by='size'
    - `max_cells`: int, default 20000000

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `part_count` and `parts`; every part has `name`, `record_count`, `variant_count`, `bytes`, `sha256` and its content (`text`, or `data_base64` for binary formats); with `options['output_dir']` each part is written to `part_0001_<name><ext>`, `part_0002_...` and `output_path` is returned instead of the content. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import split_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = split_ped(PED, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "part_count": 1,
    "parts": [
      {
        "name": "F1",
        "record_count": 3,
        "variant_count": 2,
        "bytes": 69,
        "sha256": "5083df15f1e129b12dfd6258db03b28858011c33cce5a983955ed2f3ace7a2b7",
        "text": "F1\tI1\t0\t0\t1\t2\tA\tA\tC\tG\nF1\tI2\t0\t0\t2\t1\tA\tG\t0\t0\nF1\tI3\tI1\tI2\t1\t-9\tG\tG\tC\tC\n"
      }
    ]
  }
}
```

## 91) split_phylip

- **Purpose:** Split a PHYLIP file by site blocks, taxon groups or named site partitions.

- **Method:** Give exactly one of `site_size` (consecutive blocks of sites, parts named `sites_<first>-<last>`), `taxa_size` (consecutive groups of taxa, `taxa_0001`, ...) or `partitions` (label to `start-end[,start-end...]`, 1-based inclusive site ranges, one part per label). PHYLIP files hold no trees, so the `trees_per_part` option of the NEXUS and Newick splitters is rejected here with a `ConfigurationError` that says PHYLIP has no trees. Each part is a complete PHYLIP file.

- **API:** `dnakit.io_formats.split_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PHYLIP alignment file (header `ntax nchar`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output_dir`: directory to write the parts to as part_0001_<label><ext> (optional)
    - `overwrite`: bool, default False
    - `datatype`: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing characters (gap, ?, N for nucleotides, X for protein)
    - `input_layout`: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
    - `input_strict_names`: bool, default False: input names occupy the first 10 characters of each row (otherwise the first whitespace-delimited word)
    - `site_size`: int: consecutive blocks of this many sites
    - `taxa_size`: int: consecutive groups of this many taxa
    - `partitions`: mapping label -> 'start-end[,start-end...]' (1-based, inclusive site ranges). Give exactly one of site_size, taxa_size and partitions

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `part_count` and `parts`; every part has `name`, `record_count`, `bytes`, `sha256` and its content (`text`, or `data_base64` for binary formats); with `options['output_dir']` each part is written to `part_0001_<name><ext>`, `part_0002_...` and `output_path` is returned instead of the content. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import split_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = split_phylip(PHYLIP, options={"site_size": 4}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "part_count": 2,
    "parts": [
      {
        "name": "sites_1-4",
        "record_count": 3,
        "bytes": 53,
        "sha256": "878fecc7ee77357a027fd67f53ca9a126c676fd3fb4de77519f8d7a15695fcd2",
        "text": " 3 4\nA          ACGT\nB          ACGT\nC          ACG-\n"
      },
      {
        "name": "sites_5-8",
        "record_count": 3,
        "bytes": 53,
        "sha256": "05cbfd292180c3c0d25dc5286ea79a552121d97f8d753d01a20fb56bd2c1995d",
        "text": " 3 4\nA          ACGT\nB          ACGA\nC          ACGT\n"
      }
    ]
  }
}
```

## 92) split_twobit

- **Purpose:** Split a 2bit file into 2bit files of N sequences.

- **Method:** Consecutive groups of `sequences_per_part` sequences (default 1) are written to separate version-0 files; each part is a complete 2bit file. With `sequences_per_part=1` a part is named after its sequence, otherwise parts are named `0001`, `0002`, ....

- **API:** `dnakit.io_formats.split_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a UCSC 2bit genome file, given as a path (`str` or `Path`) or as its raw `bytes` (a `str` that is not an existing path is rejected, because text cannot hold a binary file); both byte orders are accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `output_dir`: directory to write the parts to as part_0001_<label><ext> (optional)
    - `overwrite`: bool, default False
    - `sequences_per_part`: int, default 1: sequences in each output part
    - `byte_order`: 'little' (default) or 'big': byte order of the written 2bit file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `part_count` and `parts`; every part has `name`, `record_count`, `sequence_count`, `bytes`, `sha256` and its content (`text`, or `data_base64` for binary formats); with `options['output_dir']` each part is written to `part_0001_<name><ext>`, `part_0002_...` and `output_path` is returned instead of the content. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import base64
import json

from dnakit.io_formats import split_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = split_twobit(TWOBIT, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "part_count": 2,
    "parts": [
      {
        "name": "s1",
        "record_count": 1,
        "sequence_count": 1,
        "bytes": 58,
        "sha256": "2c8722ce953b5a3db475f68d7f8f16d9f09381af10b7e1472ac129e3eed9f6cf",
        "data_base64": "QydBGgAAAAABAAAAAAAAAAJzMRcAAAAKAAAAAQAAAAQAAAACAAAAAQAAAAYAAAAEAAAAAAAAAJwJwA=="
      },
      {
        "name": "s2",
        "record_count": 1,
        "sequence_count": 1,
        "bytes": 41,
        "sha256": "99bac5b233c34dd7e0e7a506a9b1b176814d4680c658ebf52792b0681b440e44",
        "data_base64": "QydBGgAAAAABAAAAAAAAAAJzMhcAAAAIAAAAAAAAAAAAAAAAAAAA9aA="
      }
    ]
  }
}
```

## 93) split_zarr

- **Purpose:** Split a Zarr store into stores of N rows along axis 0.

- **Method:** Axis 0 of all arrays is cut into consecutive blocks of `rows_per_part` rows (required) and one store per block is written to `output_dir` (required) as `part_0001_rows_<first>-<last>.zarr`, ...; all arrays with at least one dimension must have the same length on axis 0 (otherwise `ConfigurationError`; use `filter_zarr` to select compatible arrays) and arrays without a dimension are copied into every part. Splitting more than `max_elements` elements is refused.

- **API:** `dnakit.io_formats.split_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the path (`str` or `Path`) of a Zarr directory store (zarr v2 or v3; needs the optional `zarr>=3` and `numpy` packages); text content, bytes and gzip are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `output_dir`: directory that receives part_0001_<label>.zarr stores (required)
    - `overwrite`: bool, default False
    - `rows_per_part`: int (required): rows along axis 0 per part
    - `max_elements`: int, default 5000000: refuse to materialise more array elements than this

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `part_count` and `parts` (per part `name` such as `rows_0-3`, `record_count` (arrays), `row_count` and `output_path`; the parts are always written to `output_dir`). If `zarr` (>= 3) or `numpy` cannot be imported the function returns `status: "not_implemented"` with `required_dependency: "pip install 'zarr>=3' numpy"` and a `message` (`zarr and numpy are required: <import error>`, or `zarr <version> is installed but this adapter needs zarr >= 3.`) instead of a `result`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import split_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = split_zarr("demo.zarr", options={"rows_per_part": 3, "output_dir": "parts"}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "part_count": 2,
    "parts": [
      {
        "name": "rows_0-3",
        "record_count": 2,
        "row_count": 3,
        "output_path": "parts/part_0001_rows_0-3.zarr"
      },
      {
        "name": "rows_3-4",
        "record_count": 2,
        "row_count": 1,
        "output_path": "parts/part_0002_rows_3-4.zarr"
      }
    ]
  }
}
```

## 94) stream

- **Purpose:** Return one page of records from a file of any supported format, with a cursor for the next page.

- **Method:** The format is taken from `options['format']` (`map`, `newick`, `nexus`, `paf`, `ped`, `phylip`, `twobit` or `zarr`; `fasta` only as the source of `convert`) or detected, first from the path (suffixes `.map`, `.ped`, `.paf`, `.nwk`/`.newick`/`.tre`/`.tree`, `.nex`/`.nexus`/`.nxs`, `.phy`/`.phylip`, `.2bit`, `.fa`/`.fasta`/`.fna`, with a trailing `.gz` ignored; a directory is taken as zarr) and otherwise from the first 64 KiB of content (2bit signature, `#NEXUS`, a leading `(` for Newick, a leading `>` for FASTA, a first line with at least 12 tab-separated columns whose columns 2 to 4 are integers for PAF, an `ntax nchar` header for PHYLIP, four columns ending in an integer for MAP, an even number of at least six columns for PED). The call is then executed by the same engine as the format-specific function (for example `stream_map`), so the options and the `result` are those of that function. If the format cannot be detected `InputFormatError` (`IO_FORMATS_UNKNOWN_FORMAT`) is raised. A FASTA input with any operation other than `convert` raises `ConfigurationError`.

- **API:** `dnakit.io_formats.stream(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a file of any supported format, as accepted by its format-specific function (path, text or bytes). `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Besides the keys below, the options of the format-specific function of the detected format are accepted. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `start`: int, default 0: index of the first record of the page
    - `limit`: int, default 1000: records per page (max 100000)
    - `format`: name of the input format (map, newick, nexus, paf, ped, phylip, twobit, zarr); detected from the suffix/content when omitted

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: exactly the `result` of the format-specific function that was used; its `format` field tells which one (for FASTA input to `convert`, `source_format` is `fasta`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import stream

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = stream(MAP, options={"limit": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "start": 0,
    "limit": 2,
    "returned": 2,
    "total": 4,
    "next_start": 2,
    "done": false,
    "records": [
      {
        "chrom": "1",
        "id": "rs1",
        "cm": 0.0,
        "pos": 100
      },
      {
        "chrom": "1",
        "id": "rs2",
        "cm": 0.5,
        "pos": 50
      }
    ]
  }
}
```

## 95) stream_map

- **Purpose:** Return one page of variants from a PLINK `.map` file, with a cursor for the next page.

- **Method:** The file is parsed (bounded by `max_bytes`) and the records `[start, start + limit)` are returned; call again with `start=next_start` until `done` is true.

- **API:** `dnakit.io_formats.stream_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.map` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `start`: int, default 0: index of the first record of the page
    - `limit`: int, default 1000: records per page (max 100000)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `start`, `limit`, `returned`, `total`, `next_start`, `done` and `records` (each record is `{chrom, id, cm, pos}`); `next_start` is the value to pass as `options['start']` for the next page and is `null` on the last page, where `done` is `true`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import stream_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = stream_map(MAP, options={"start": 1, "limit": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "start": 1,
    "limit": 2,
    "returned": 2,
    "total": 4,
    "next_start": 3,
    "done": false,
    "records": [
      {
        "chrom": "1",
        "id": "rs2",
        "cm": 0.5,
        "pos": 50
      },
      {
        "chrom": "2",
        "id": "rs3",
        "cm": 0.0,
        "pos": 10
      }
    ]
  }
}
```

## 96) stream_newick

- **Purpose:** Return one page of trees from a Newick file, with a cursor for the next page.

- **Method:** The file is parsed (bounded by `max_bytes` and `max_nodes`) and the trees `[start, start + limit)` are returned; call again with `start=next_start` until `done` is true.

- **API:** `dnakit.io_formats.stream_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a Newick tree file (one or more `;`-terminated trees), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `start`: int, default 0: index of the first record of the page
    - `limit`: int, default 1000: records per page (max 100000)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `start`, `limit`, `returned`, `total`, `next_start`, `done` and `records` (each record is `{index, newick, tips, nodes}` with the normalised Newick string, the tip names and the flat node table `{id, parent, name, length}` in preorder); `next_start` is the value to pass as `options['start']` for the next page and is `null` on the last page, where `done` is `true`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import stream_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = stream_newick(TREES, options={"start": 1}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "start": 1,
    "limit": 1000,
    "returned": 1,
    "total": 2,
    "next_start": null,
    "done": true,
    "records": [
      {
        "index": 1,
        "newick": "(A,B,(C,D));",
        "tips": [
          "A",
          "B",
          "C",
          "D"
        ],
        "nodes": [
          {
            "id": 0,
            "parent": null,
            "name": null,
            "length": null
          },
          {
            "id": 1,
            "parent": 0,
            "name": "A",
            "length": null
          },
          {
            "id": 2,
            "parent": 0,
            "name": "B",
            "length": null
          },
          {
            "id": 3,
            "parent": 0,
            "name": null,
            "length": null
          },
          {
            "id": 4,
            "parent": 3,
            "name": "C",
            "length": null
          },
          {
            "id": 5,
            "parent": 3,
            "name": "D",
            "length": null
          }
        ]
      }
    ]
  }
}
```

## 97) stream_nexus

- **Purpose:** Return one page of sequences (then trees) from a NEXUS file, with a cursor for the next page.

- **Method:** The file is parsed (bounded by `max_bytes`) and the records `[start, start + limit)` of the combined list (matrix rows first, then trees) are returned; call again with `start=next_start` until `done` is true.

- **API:** `dnakit.io_formats.stream_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a NEXUS file (`#NEXUS`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `start`: int, default 0: index of the first record of the page
    - `limit`: int, default 1000: records per page (max 100000)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `start`, `limit`, `returned`, `total`, `next_start`, `done` and `records` (sequence records `{kind: "sequence", name, sequence}` followed by tree records `{kind: "tree", name, rooted, newick}`); `next_start` is the value to pass as `options['start']` for the next page and is `null` on the last page, where `done` is `true`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import stream_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = stream_nexus(NEXUS, options={"limit": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "start": 0,
    "limit": 2,
    "returned": 2,
    "total": 4,
    "next_start": 2,
    "done": false,
    "records": [
      {
        "kind": "sequence",
        "name": "A",
        "sequence": "ACGTACGT"
      },
      {
        "kind": "sequence",
        "name": "B",
        "sequence": "ACGTACGA"
      }
    ]
  }
}
```

## 98) stream_paf

- **Purpose:** Return one page of alignments from a PAF file, with a cursor for the next page.

- **Method:** The file is parsed (bounded by `max_bytes`) and the records `[start, start + limit)` are returned; call again with `start=next_start` until `done` is true.

- **API:** `dnakit.io_formats.stream_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PAF (Pairwise mApping Format) file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `start`: int, default 0: index of the first record of the page
    - `limit`: int, default 1000: records per page (max 100000)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `start`, `limit`, `returned`, `total`, `next_start`, `done` and `records` (each record has the 12 PAF fields plus `tags` (tag to `{type, value}`)); `next_start` is the value to pass as `options['start']` for the next page and is `null` on the last page, where `done` is `true`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import stream_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = stream_paf(PAF, options={"limit": 1}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "start": 0,
    "limit": 1,
    "returned": 1,
    "total": 2,
    "next_start": 1,
    "done": false,
    "records": [
      {
        "query_name": "q1",
        "query_length": 100,
        "query_start": 0,
        "query_end": 90,
        "strand": "+",
        "target_name": "t1",
        "target_length": 1000,
        "target_start": 100,
        "target_end": 190,
        "matches": 80,
        "block_length": 90,
        "mapq": 60,
        "tags": {
          "tp": {
            "type": "A",
            "value": "P"
          },
          "cg": {
            "type": "Z",
            "value": "90M"
          }
        }
      }
    ]
  }
}
```

## 99) stream_ped

- **Purpose:** Return one page of individuals from a PLINK `.ped` file, with a cursor for the next page.

- **Method:** The file is parsed (bounded by `max_bytes` and `max_cells`) and the records `[start, start + limit)` are returned; call again with `start=next_start` until `done` is true.

- **API:** `dnakit.io_formats.stream_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.ped` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `start`: int, default 0: index of the first record of the page
    - `limit`: int, default 1000: records per page (max 100000)
    - `max_cells`: int, default 20000000

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `start`, `limit`, `returned`, `total`, `next_start`, `done` and `records` (each record is `{fid, iid, pat, mat, sex, pheno, genotypes}` with `genotypes` a list of `[allele1, allele2]` pairs); `next_start` is the value to pass as `options['start']` for the next page and is `null` on the last page, where `done` is `true`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import stream_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = stream_ped(PED, options={"limit": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "start": 0,
    "limit": 2,
    "returned": 2,
    "total": 3,
    "next_start": 2,
    "done": false,
    "records": [
      {
        "fid": "F1",
        "iid": "I1",
        "pat": "0",
        "mat": "0",
        "sex": 1,
        "pheno": "2",
        "genotypes": [
          [
            "A",
            "A"
          ],
          [
            "C",
            "G"
          ]
        ]
      },
      {
        "fid": "F1",
        "iid": "I2",
        "pat": "0",
        "mat": "0",
        "sex": 2,
        "pheno": "1",
        "genotypes": [
          [
            "A",
            "G"
          ],
          [
            "0",
            "0"
          ]
        ]
      }
    ]
  }
}
```

## 100) stream_phylip

- **Purpose:** Return one page of sequences from a PHYLIP file, with a cursor for the next page.

- **Method:** The file is parsed (bounded by `max_bytes`) and the sequences `[start, start + limit)` are returned; call again with `start=next_start` until `done` is true.

- **API:** `dnakit.io_formats.stream_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PHYLIP alignment file (header `ntax nchar`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `start`: int, default 0: index of the first record of the page
    - `limit`: int, default 1000: records per page (max 100000)
    - `datatype`: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing characters (gap, ?, N for nucleotides, X for protein)
    - `input_layout`: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
    - `input_strict_names`: bool, default False: input names occupy the first 10 characters of each row (otherwise the first whitespace-delimited word)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `start`, `limit`, `returned`, `total`, `next_start`, `done` and `records` (each record is `{name, sequence}`); `next_start` is the value to pass as `options['start']` for the next page and is `null` on the last page, where `done` is `true`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import stream_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = stream_phylip(PHYLIP, options={"start": 1, "limit": 1}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "start": 1,
    "limit": 1,
    "returned": 1,
    "total": 3,
    "next_start": 2,
    "done": false,
    "records": [
      {
        "name": "B",
        "sequence": "ACGTACGA"
      }
    ]
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

- **Purpose:** Return one page of sequence windows from a 2bit file, with a cursor for the next page.

- **Method:** Every sequence is cut into windows of `window` bases (default 10000) and the windows of all sequences, in file order, form the record list; `start` and `limit` count windows, and only the windows of the page are decoded (at most `max_bases` bases per page).

- **API:** `dnakit.io_formats.stream_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a UCSC 2bit genome file, given as a path (`str` or `Path`) or as its raw `bytes` (a `str` that is not an existing path is rejected, because text cannot hold a binary file); both byte orders are accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `start`: int, default 0: index of the first record of the page
    - `limit`: int, default 1000: records per page (max 100000)
    - `window`: int, default 10000: bases per window (the unit that is paged)
    - `soft_mask`: bool, default True
    - `max_bases`: int, default 50000000: bases per page

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `start`, `limit`, `returned`, `total` (number of windows), `window`, `next_start`, `done` and `records`; each record is `{name, start, end, sequence}` for one window (`next_start` is `null` on the last page). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import base64
import json

from dnakit.io_formats import stream_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = stream_twobit(TWOBIT, options={"window": 6, "limit": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "start": 0,
    "limit": 2,
    "returned": 2,
    "total": 4,
    "window": 6,
    "next_start": 2,
    "done": false,
    "records": [
      {
        "name": "s1",
        "start": 0,
        "end": 6,
        "sequence": "ACGTNN"
      },
      {
        "name": "s1",
        "start": 6,
        "end": 10,
        "sequence": "acgt"
      }
    ]
  }
}
```

## 103) stream_zarr

- **Purpose:** Return one page of rows from one array of a Zarr store, with a cursor for the next page.

- **Method:** Pages through the rows along axis 0 of one array: each record is a row as a nested list and `start`/`limit` count rows. `array` is required when the store has several arrays. At most `max_elements` elements are read per page.

- **API:** `dnakit.io_formats.stream_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the path (`str` or `Path`) of a Zarr directory store (zarr v2 or v3; needs the optional `zarr>=3` and `numpy` packages); text content, bytes and gzip are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `start`: int, default 0: index of the first record of the page
    - `limit`: int, default 1000: records per page (max 100000)
    - `array`: array path to page through (required when the store has several arrays)
    - `max_elements`: int, default 5000000: maximum elements per page

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `array`, `start`, `limit`, `returned`, `total`, `next_start`, `done` and `records` (rows of the array as nested lists); `next_start` is the value to pass as `options['start']` for the next page and is `null` on the last page, where `done` is `true`. If `zarr` (>= 3) or `numpy` cannot be imported the function returns `status: "not_implemented"` with `required_dependency: "pip install 'zarr>=3' numpy"` and a `message` (`zarr and numpy are required: <import error>`, or `zarr <version> is installed but this adapter needs zarr >= 3.`) instead of a `result`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import stream_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = stream_zarr("demo.zarr", options={"array": "dosage", "start": 2, "limit": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "array": "dosage",
    "start": 2,
    "limit": 2,
    "returned": 2,
    "total": 4,
    "next_start": null,
    "done": true,
    "records": [
      [
        0,
        0
      ],
      [
        2,
        1
      ]
    ]
  }
}
```

## 104) summarize

- **Purpose:** Summarise a file of any supported format.

- **Method:** The format is taken from `options['format']` (`map`, `newick`, `nexus`, `paf`, `ped`, `phylip`, `twobit` or `zarr`; `fasta` only as the source of `convert`) or detected, first from the path (suffixes `.map`, `.ped`, `.paf`, `.nwk`/`.newick`/`.tre`/`.tree`, `.nex`/`.nexus`/`.nxs`, `.phy`/`.phylip`, `.2bit`, `.fa`/`.fasta`/`.fna`, with a trailing `.gz` ignored; a directory is taken as zarr) and otherwise from the first 64 KiB of content (2bit signature, `#NEXUS`, a leading `(` for Newick, a leading `>` for FASTA, a first line with at least 12 tab-separated columns whose columns 2 to 4 are integers for PAF, an `ntax nchar` header for PHYLIP, four columns ending in an integer for MAP, an even number of at least six columns for PED). The call is then executed by the same engine as the format-specific function (for example `summarize_twobit`), so the options and the `result` are those of that function. If the format cannot be detected `InputFormatError` (`IO_FORMATS_UNKNOWN_FORMAT`) is raised. A FASTA input with any operation other than `convert` raises `ConfigurationError`.

- **API:** `dnakit.io_formats.summarize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a file of any supported format, as accepted by its format-specific function (path, text or bytes). `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Besides the keys below, the options of the format-specific function of the detected format are accepted. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `format`: name of the input format (map, newick, nexus, paf, ped, phylip, twobit, zarr); detected from the suffix/content when omitted

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: exactly the `result` of the format-specific function that was used; its `format` field tells which one (for FASTA input to `convert`, `source_format` is `fasta`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import base64
import json

from dnakit.io_formats import summarize, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = summarize(TWOBIT, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "sequence_count": 2,
    "total_length": 18,
    "n_bases": 2,
    "soft_masked_bases": 4,
    "gc_fraction": 0.5,
    "version": 0,
    "byte_order": "little",
    "sequences": [
      {
        "name": "s1",
        "length": 10,
        "A": 2,
        "C": 2,
        "G": 2,
        "T": 2,
        "N": 2,
        "soft_masked": 4,
        "gc_fraction": 0.5
      },
      {
        "name": "s2",
        "length": 8,
        "A": 2,
        "C": 2,
        "G": 2,
        "T": 2,
        "N": 0,
        "soft_masked": 0,
        "gc_fraction": 0.5
      }
    ]
  }
}
```

## 105) summarize_map

- **Purpose:** Summarise a PLINK `.map` file: variant and chromosome counts, bp ranges and duplicate IDs.

- **Method:** Counts variants per chromosome with their minimum and maximum bp (chromosomes listed in natural order 1..22, X, Y, XY, M/MT, others), the number of IDs that occur more than once, the range of cM values, whether all cM values are 0, and whether the file is already sorted by chromosome and position.

- **API:** `dnakit.io_formats.summarize_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.map` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `variant_count`, `chromosome_count`, `chromosomes` (per chromosome `variant_count`, `min_pos`, `max_pos`), `duplicate_id_count`, `cm_min`, `cm_max`, `all_cm_zero` and `sorted`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import summarize_map

MAP = (
    "1\trs1\t0\t100\n"
    "1\trs2\t0.5\t50\n"
    "2\trs3\t0\t10\n"
    "X\trs4\t0\t5\n"
)

envelope = summarize_map(MAP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "variant_count": 4,
    "chromosome_count": 3,
    "chromosomes": {
      "1": {
        "variant_count": 2,
        "min_pos": 50,
        "max_pos": 100
      },
      "2": {
        "variant_count": 1,
        "min_pos": 10,
        "max_pos": 10
      },
      "X": {
        "variant_count": 1,
        "min_pos": 5,
        "max_pos": 5
      }
    },
    "duplicate_id_count": 0,
    "cm_min": 0.0,
    "cm_max": 0.5,
    "all_cm_zero": false,
    "sorted": false
  }
}
```

## 106) summarize_newick

- **Purpose:** Summarise a Newick file: per-tree node counts, depth, height, total branch length and resolution.

- **Method:** Per tree: tip, internal and total node counts; root degree; `rooted` (the root has exactly two children); `fully_resolved` flag; depth in edges; height = longest root-to-tip path length (a missing length counts as 0); total branch length (`null` when no branch has a length); counts of negative branches, unnamed tips and duplicated tip names.

- **API:** `dnakit.io_formats.summarize_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a Newick tree file (one or more `;`-terminated trees), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `tree_count`, `trees` (per tree `tip_count`, `internal_node_count`, `node_count`, `root_degree`, `rooted`, `fully_resolved`, `max_depth_edges`, `height`, `total_branch_length`, `has_branch_lengths`, `negative_branch_count`, `unnamed_tip_count`, `duplicate_tip_count`), `tip_count_min` and `tip_count_max`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import summarize_newick

TREES = "((A:1,B:2):0.5,(C:1,D:1):0.5);\n(A,B,(C,D));\n"

envelope = summarize_newick(TREES, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "tree_count": 2,
    "trees": [
      {
        "tip_count": 4,
        "internal_node_count": 3,
        "node_count": 7,
        "root_degree": 2,
        "rooted": true,
        "fully_resolved": true,
        "max_depth_edges": 2,
        "height": 2.5,
        "total_branch_length": 6.0,
        "has_branch_lengths": true,
        "negative_branch_count": 0,
        "unnamed_tip_count": 0,
        "duplicate_tip_count": 0
      },
      {
        "tip_count": 4,
        "internal_node_count": 2,
        "node_count": 6,
        "root_degree": 3,
        "rooted": false,
        "fully_resolved": true,
        "max_depth_edges": 2,
        "height": 0.0,
        "total_branch_length": null,
        "has_branch_lengths": false,
        "negative_branch_count": 0,
        "unnamed_tip_count": 0,
        "duplicate_tip_count": 0
      }
    ],
    "tip_count_min": 4,
    "tip_count_max": 4
  }
}
```

## 107) summarize_nexus

- **Purpose:** Summarise a NEXUS file: alignment statistics (constant, variable and parsimony-informative sites, composition) and tree statistics.

- **Method:** Alignment statistics as for `summarize_phylip` (a column is constant when its unambiguous characters show one state, variable when they show several, parsimony-informative when at least two states each occur in at least two taxa; gaps, missing symbols and ambiguity codes are ignored, `U` is counted as `T` for nucleotides), plus the `summarize_newick` statistics of each tree.

- **API:** `dnakit.io_formats.summarize_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a NEXUS file (`#NEXUS`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `taxon_count`, `site_count`, `datatype`, `equal_length`, `constant_sites`, `variable_sites`, `parsimony_informative_sites`, `all_missing_sites`, `total_characters`, `missing_fraction`, `taxon_missing_fraction` (`min`, `mean`, `max`), `duplicate_sequence_count`, `composition`, `gc_content` (nucleotides) and, when trees exist, `tree_count` and `trees`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import summarize_nexus

NEXUS = """#NEXUS
BEGIN TAXA;
  DIMENSIONS NTAX=3;
  TAXLABELS A B 'C d';
END;
BEGIN CHARACTERS;
  DIMENSIONS NCHAR=8;
  FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=NO;
  MATRIX
    A ACGT ACGT
    B ACGTACGA
    'C d' ACG-ACGT
  ;
END;
BEGIN TREES;
  TRANSLATE 1 A, 2 B, 3 'C d';
  TREE t1 = [&R] ((1:1,2:1):1,3:2);
END;
"""

envelope = summarize_nexus(NEXUS, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "taxon_count": 3,
    "site_count": 8,
    "datatype": "dna",
    "equal_length": true,
    "constant_sites": 7,
    "variable_sites": 1,
    "parsimony_informative_sites": 0,
    "all_missing_sites": 0,
    "total_characters": 24,
    "missing_fraction": 0.041666666666666664,
    "taxon_missing_fraction": {
      "min": 0.0,
      "mean": 0.041666666666666664,
      "max": 0.125
    },
    "duplicate_sequence_count": 0,
    "composition": {
      "-": 1,
      "A": 7,
      "C": 6,
      "G": 6,
      "T": 4
    },
    "gc_content": 0.5217391304347826,
    "tree_count": 1,
    "trees": [
      {
        "tip_count": 3,
        "internal_node_count": 2,
        "node_count": 5,
        "root_degree": 2,
        "rooted": true,
        "fully_resolved": true,
        "max_depth_edges": 2,
        "height": 2.0,
        "total_branch_length": 5.0,
        "has_branch_lengths": true,
        "negative_branch_count": 0,
        "unnamed_tip_count": 0,
        "duplicate_tip_count": 0
      }
    ]
  }
}
```

## 108) summarize_paf

- **Purpose:** Summarise a PAF file: alignment counts, identity, aligned bases, strands and the mapping-quality histogram.

- **Method:** Counts alignments and distinct queries/targets; aligned bases = sum of block lengths; per-alignment identity = matches / block_length, reported as mean, minimum and maximum, and `overall_identity` = total matches / total block length (base-weighted); `mapq` is the count of alignments per mapping-quality value. Empty input gives `null` for the ratios.

- **API:** `dnakit.io_formats.summarize_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PAF (Pairwise mApping Format) file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `alignment_count`, `query_count`, `target_count`, `strand` (`+`, `-`), `aligned_bases`, `matching_bases`, `overall_identity`, `mean_identity`, `min_identity`, `max_identity`, `mapq` (value to count) and `mean_mapq`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import summarize_paf

PAF = (
    "q1\t100\t0\t90\t+\tt1\t1000\t100\t190\t80\t90\t60\ttp:A:P\tcg:Z:90M\n"
    "q2\t50\t0\t50\t-\tt2\t500\t0\t50\t40\t50\t10\n"
)

envelope = summarize_paf(PAF, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "alignment_count": 2,
    "query_count": 2,
    "target_count": 2,
    "strand": {
      "+": 1,
      "-": 1
    },
    "aligned_bases": 140,
    "matching_bases": 120,
    "overall_identity": 0.8571428571428571,
    "mean_identity": 0.8444444444444444,
    "min_identity": 0.8,
    "max_identity": 0.8888888888888888,
    "mapq": {
      "10": 1,
      "60": 1
    },
    "mean_mapq": 35.0
  }
}
```

## 109) summarize_ped

- **Purpose:** Summarise a PLINK `.ped` file: individual, family and variant counts, sex, phenotype classes and missing-genotype rates.

- **Method:** Counts individuals, families, variants and founders (father and mother both `0`); sex is classified as male (1), female (2) or unknown (anything else); phenotype as case (2), control (1), missing (0 or -9) or quantitative (any other value). A genotype call counts as missing when either allele is `0`; `genotype_missing_rate` is missing calls / all calls and `individual_missing_rate` gives the minimum, mean and maximum per-individual rate.

- **API:** `dnakit.io_formats.summarize_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.ped` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `max_cells`: int, default 20000000

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `individual_count`, `family_count`, `variant_count`, `founder_count`, `sex` (`male`, `female`, `unknown`), `phenotype` (`case`, `control`, `missing`, `quantitative`), `genotype_missing_rate` and `individual_missing_rate` (`min`, `mean`, `max`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import summarize_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 2 1 A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = summarize_ped(PED, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "individual_count": 3,
    "family_count": 1,
    "variant_count": 2,
    "founder_count": 2,
    "sex": {
      "male": 2,
      "female": 1,
      "unknown": 0
    },
    "phenotype": {
      "case": 1,
      "control": 1,
      "missing": 1,
      "quantitative": 0
    },
    "genotype_missing_rate": 0.16666666666666666,
    "individual_missing_rate": {
      "min": 0.0,
      "mean": 0.16666666666666666,
      "max": 0.5
    }
  }
}
```

## 110) summarize_phylip

- **Purpose:** Summarise a PHYLIP file: alignment statistics (constant, variable and parsimony-informative sites, composition, missing data).

- **Method:** A column is constant when its unambiguous characters show one state, variable when they show several, and parsimony-informative when at least two different unambiguous states each occur in at least two taxa (gaps, missing symbols and ambiguity codes are ignored; for DNA/RNA `U` is counted as `T`; columns without any unambiguous character are `all_missing_sites`). `gc_content` = (G + C) / (A + C + G + T + U) over the sequence characters, for nucleotide data; `missing_fraction` and the per-taxon fractions count gap, missing and ambiguity characters of the `datatype`; duplicate sequences are counted.

- **API:** `dnakit.io_formats.summarize_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PHYLIP alignment file (header `ntax nchar`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456
    - `datatype`: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing characters (gap, ?, N for nucleotides, X for protein)
    - `input_layout`: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
    - `input_strict_names`: bool, default False: input names occupy the first 10 characters of each row (otherwise the first whitespace-delimited word)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `taxon_count`, `site_count`, `datatype`, `equal_length`, `constant_sites`, `variable_sites`, `parsimony_informative_sites`, `all_missing_sites`, `total_characters`, `missing_fraction`, `taxon_missing_fraction` (`min`, `mean`, `max`), `duplicate_sequence_count`, `composition` (character to count) and `gc_content` (nucleotides). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import summarize_phylip

PHYLIP = " 3 8\nA  ACGTACGT\nB  ACGTACGA\nC  ACG-ACGT\n"

envelope = summarize_phylip(PHYLIP, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "taxon_count": 3,
    "site_count": 8,
    "datatype": "dna",
    "equal_length": true,
    "constant_sites": 7,
    "variable_sites": 1,
    "parsimony_informative_sites": 0,
    "all_missing_sites": 0,
    "total_characters": 24,
    "missing_fraction": 0.041666666666666664,
    "taxon_missing_fraction": {
      "min": 0.0,
      "mean": 0.041666666666666664,
      "max": 0.125
    },
    "duplicate_sequence_count": 0,
    "composition": {
      "-": 1,
      "A": 7,
      "C": 6,
      "G": 6,
      "T": 4
    },
    "gc_content": 0.5217391304347826
  }
}
```

## 111) summarize_twobit

- **Purpose:** Summarise a 2bit file: per-sequence length, base composition, N and soft-masked bases, GC fraction.

- **Method:** Per sequence: length, counts of A, C, G, T (case-insensitive) and N, the number of soft-masked (lower-case) bases and `gc_fraction = (G + C) / (A + C + G + T)`; the file totals add the sequence lengths, N bases and masked bases and compute `gc_fraction` over all sequences.

- **API:** `dnakit.io_formats.summarize_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a UCSC 2bit genome file, given as a path (`str` or `Path`) or as its raw `bytes` (a `str` that is not an existing path is rejected, because text cannot hold a binary file); both byte orders are accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):


- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `sequence_count`, `total_length`, `n_bases`, `soft_masked_bases`, `gc_fraction`, `version`, `byte_order` and `sequences` (per sequence `name`, `length`, `A`, `C`, `G`, `T`, `N`, `soft_masked`, `gc_fraction`). Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import base64
import json

from dnakit.io_formats import summarize_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = summarize_twobit(TWOBIT, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "sequence_count": 2,
    "total_length": 18,
    "n_bases": 2,
    "soft_masked_bases": 4,
    "gc_fraction": 0.5,
    "version": 0,
    "byte_order": "little",
    "sequences": [
      {
        "name": "s1",
        "length": 10,
        "A": 2,
        "C": 2,
        "G": 2,
        "T": 2,
        "N": 2,
        "soft_masked": 4,
        "gc_fraction": 0.5
      },
      {
        "name": "s2",
        "length": 8,
        "A": 2,
        "C": 2,
        "G": 2,
        "T": 2,
        "N": 0,
        "soft_masked": 0,
        "gc_fraction": 0.5
      }
    ]
  }
}
```

## 112) summarize_zarr

- **Purpose:** Summarise a Zarr store: array shapes, dtypes, chunking and totals.

- **Method:** Counts groups and arrays and totals the elements, uncompressed bytes and initialised chunks; one entry per array with its metadata. No chunk data is read.

- **API:** `dnakit.io_formats.summarize_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the path (`str` or `Path`) of a Zarr directory store (zarr v2 or v3; needs the optional `zarr>=3` and `numpy` packages); text content, bytes and gzip are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):


- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `path`, `kind` (`"group"` or `"array"`), `array_count`, `group_count`, `total_elements`, `total_nbytes`, `initialized_chunks`, `total_chunks`, `attrs` and `arrays` (per array shape, dtype, chunks, chunk grid, chunk counts, compressors, zarr format, `nbytes`, `attrs`). If `zarr` (>= 3) or `numpy` cannot be imported the function returns `status: "not_implemented"` with `required_dependency: "pip install 'zarr>=3' numpy"` and a `message` (`zarr and numpy are required: <import error>`, or `zarr <version> is installed but this adapter needs zarr >= 3.`) instead of a `result`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import summarize_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = summarize_zarr("demo.zarr", verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "path": "demo.zarr",
    "kind": "group",
    "array_count": 2,
    "group_count": 0,
    "total_elements": 12,
    "total_nbytes": 96,
    "initialized_chunks": 2,
    "total_chunks": 2,
    "attrs": {},
    "arrays": {
      "dosage": {
        "shape": [
          4,
          2
        ],
        "dtype": "int64",
        "chunks": [
          4,
          2
        ],
        "chunk_grid": [
          1,
          1
        ],
        "nchunks": 1,
        "nchunks_initialized": 1,
        "compressors": [
          "ZstdCodec"
        ],
        "zarr_format": 3,
        "nbytes": 64,
        "attrs": {}
      },
      "pos": {
        "shape": [
          4
        ],
        "dtype": "int64",
        "chunks": [
          4
        ],
        "chunk_grid": [
          1
        ],
        "nchunks": 1,
        "nchunks_initialized": 1,
        "compressors": [
          "ZstdCodec"
        ],
        "zarr_format": 3,
        "nbytes": 32,
        "attrs": {}
      }
    }
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

- **Purpose:** Validate a PLINK `.map` file and report its errors and warnings.

- **Method:** Each non-blank line needs four whitespace-separated columns: chromosome, variant ID, genetic distance in cM (finite number) and base-pair position (integer >= 0). Every problem of a line is reported as its own entry `line N: field 'NAME' (column C) has invalid value 'VALUE': WHY` (for example a cM of `x` gives `field 'cM' (column 3) has invalid value 'x': not a finite number`, and a line with a bad cM and a bad bp gives two entries); a line with the wrong number of columns gives one entry `line N: expected 4 columns ...`. Duplicate variant IDs are errors; duplicate (chromosome, bp) pairs are warnings.

- **API:** `dnakit.io_formats.validate_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.map` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `valid` (true when `errors` is empty), `errors`, `warnings` (lists of `line N: ...` messages) and `variant_count` (variants that parsed). Content problems are reported in `errors`/`warnings`, not raised; an input that cannot be read at all (missing file, invalid gzip or UTF-8) still raises `InputFormatError`; invalid options raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import validate_map

BAD = "1\trs1\t0\t100\n1\trs1\t0\t200\n2\trs3\tx\tabc\n2\trs9\t0\t10\n"

envelope = validate_map(BAD, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "valid": false,
    "errors": [
      "line 3: field 'cM' (column 3) has invalid value 'x': not a finite number",
      "line 3: field 'bp' (column 4) has invalid value 'abc': not an integer"
    ],
    "warnings": []
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

- **Purpose:** Validate a Newick tree file and report its errors and warnings.

- **Method:** Each `;`-terminated tree must have balanced parentheses, well-formed labels and finite numeric branch lengths (a syntax error is reported as `tree N: ...`). Duplicate tip names and tips without a name are errors; negative branch lengths and single-child (unary) internal nodes are warnings.

- **API:** `dnakit.io_formats.validate_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a Newick tree file (one or more `;`-terminated trees), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `valid`, `errors`, `warnings` (lists of `tree N: ...` messages), `tree_count` and `tip_count` (tips over all trees). Content problems are reported in `errors`/`warnings`, not raised; an input that cannot be read at all (missing file, invalid gzip or UTF-8) still raises `InputFormatError`; invalid options raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import validate_newick

BAD = "((A:1,A:2):0.5,(C:-1,D:1));\n"

envelope = validate_newick(BAD, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "valid": false,
    "errors": [
      "tree 1: 1 duplicated tip name(s)"
    ],
    "warnings": [
      "tree 1: 1 negative branch length(s)"
    ],
    "tree_count": 1,
    "tip_count": 4
  }
}
```

## 117) validate_nexus

- **Purpose:** Validate a NEXUS file and report its errors and warnings.

- **Method:** Checks the `#NEXUS` header, the block structure, the declared NTAX/NCHAR against the matrix, equal row lengths, unique taxa and the alphabet of nucleotide/protein data; tree tips that are absent from the matrix are reported as warnings.

- **API:** `dnakit.io_formats.validate_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a NEXUS file (`#NEXUS`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `valid`, `errors`, `warnings`, `taxon_count`, `site_count` and `tree_count`. Content problems are reported in `errors`/`warnings`, not raised; an input that cannot be read at all (missing file, invalid gzip or UTF-8) still raises `InputFormatError`; invalid options raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import validate_nexus

BAD = "#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=2 NCHAR=5;\n  FORMAT DATATYPE=DNA;\n  MATRIX\n    A ACGTA\n    B ACG\n  ;\nEND;\n"

envelope = validate_nexus(BAD, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "valid": false,
    "errors": [
      "sequences have different lengths [3, 5]; an alignment needs equal lengths",
      "header declares 5 characters but sequences have length [3, 5]"
    ],
    "warnings": [],
    "taxon_count": 2,
    "site_count": 5,
    "tree_count": 0
  }
}
```

## 118) validate_paf

- **Purpose:** Validate a PAF file and report its errors and warnings.

- **Method:** Each non-blank line needs the 12 mandatory tab-separated columns with integer lengths/coordinates/counts, `start <= end <= length` on both query and target, `matches <= block_length`, strand `+` or `-` and `mapq` in [0, 255]; optional columns must look like `TAG:TYPE:VALUE`. Every problem of a line is reported as its own entry (`line N: field 'NAME' (column C) has invalid value 'VALUE': WHY` for a single field, or `line N: query coordinates must satisfy ...` for a relation between fields), so one line can produce several entries; a line with fewer than 12 columns gives one entry.

- **API:** `dnakit.io_formats.validate_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PAF (Pairwise mApping Format) file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `valid`, `errors`, `warnings` (lists of `line N: ...` messages) and, when the file is valid, `alignment_count`. Content problems are reported in `errors`/`warnings`, not raised; an input that cannot be read at all (missing file, invalid gzip or UTF-8) still raises `InputFormatError`; invalid options raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import validate_paf

BAD = "q1\t100\t0\t90\t*\tt1\t1000\t100\t190\t95\t90\t300\n"

envelope = validate_paf(BAD, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "valid": false,
    "errors": [
      "line 1: field 'strand' (column 5) has invalid value '*': must be '+' or '-'",
      "line 1: matches must satisfy 0 <= matches <= block_length",
      "line 1: field 'mapq' (column 12) has invalid value '300': must be in [0, 255]"
    ],
    "warnings": []
  }
}
```

## 119) validate_ped

- **Purpose:** Validate a PLINK `.ped` file and report its errors and warnings.

- **Method:** Rows need `6 + 2k` whitespace-separated columns with the same `k` on every row (otherwise `line N: has ... variants but earlier rows have ...`), an integer sex and a numeric phenotype. Every problem of a row is reported as its own entry `line N: field 'NAME' (column C) has invalid value 'VALUE': WHY` (so a row with a bad sex and a bad phenotype gives two entries, and a wrong variant count does not hide them). Duplicate (family, individual) pairs are errors; a sex code outside {0, 1, 2} and a parent that is not an individual of the same family are warnings; allele `0` marks a missing call. With `variant_count` the file's variant count must equal the given number (otherwise an error).

- **API:** `dnakit.io_formats.validate_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PLINK `.ped` file, given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `variant_count`: int: expected number of variants (e.g. the line count of the matching .map)
    - `max_cells`: int, default 20000000

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `valid`, `errors`, `warnings` (lists of `line N: ...` messages), `individual_count` and `variant_count`. Content problems are reported in `errors`/`warnings`, not raised; an input that cannot be read at all (missing file, invalid gzip or UTF-8) still raises `InputFormatError`; invalid options raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import validate_ped

PED = (
    "F1 I1 0 0 1 2 A A C G\n"
    "F1 I2 0 0 two 1x A G 0 0\n"
    "F1 I3 I1 I2 1 -9 G G C C\n"
)

envelope = validate_ped(PED, options={"variant_count": 2}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "valid": false,
    "errors": [
      "line 2: field 'sex' (column 5) has invalid value 'two': not an integer",
      "line 2: field 'phenotype' (column 6) has invalid value '1x': not a finite number"
    ],
    "warnings": []
  }
}
```

## 120) validate_phylip

- **Purpose:** Validate a PHYLIP alignment file and report its errors and warnings.

- **Method:** Checks the header counts against the rows (if they do not form `ntax` taxa of exactly `nchar` characters in a sequential or interleaved layout this is one error), equal sequence lengths, unique taxon names and, for the chosen `datatype`, that every character belongs to the alphabet.

- **API:** `dnakit.io_formats.validate_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a PHYLIP alignment file (header `ntax nchar`), given as a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse inputs larger than this
    - `datatype`: 'dna' (default), 'rna', 'protein' or 'any': alphabet used for validation and for counting missing characters (gap, ?, N for nucleotides, X for protein)
    - `input_layout`: 'auto' (default), 'sequential' or 'interleaved': layout of the PHYLIP input
    - `input_strict_names`: bool, default False: input names occupy the first 10 characters of each row (otherwise the first whitespace-delimited word)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `valid`, `errors`, `warnings` and, when the file parses, `taxon_count`, `site_count` and `tree_count` (always 0 for PHYLIP). Content problems are reported in `errors`/`warnings`, not raised; an input that cannot be read at all (missing file, invalid gzip or UTF-8) still raises `InputFormatError`; invalid options raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import validate_phylip

BAD = " 3 4\nA  ACGT\nA  ACGA\nC  ACJT\n"

envelope = validate_phylip(BAD, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "valid": false,
    "errors": [
      "duplicate taxon names: ['A']",
      "C: characters ['J'] are not valid for datatype 'dna' (set the datatype option if this is a different alphabet)"
    ],
    "warnings": [],
    "taxon_count": 3,
    "site_count": 4,
    "tree_count": 0
  }
}
```

## 121) validate_twobit

- **Purpose:** Validate a UCSC 2bit file and report its errors and warnings.

- **Method:** Checks the signature (both byte orders), version 0, the reserved fields, the index offsets, unique non-empty names, that N and mask blocks are sorted, non-empty, non-overlapping and inside the sequence, and that the packed data fits in the file. A file that is not a 2bit file at all is reported as an error (`sequence_count` 0).

- **API:** `dnakit.io_formats.validate_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a UCSC 2bit genome file, given as a path (`str` or `Path`) or as its raw `bytes` (a `str` that is not an existing path is rejected, because text cannot hold a binary file); both byte orders are accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):


- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `valid`, `errors`, `warnings`, `sequence_count` and, for a readable header, `byte_order` and `file_size`. Content problems are reported in `errors`/`warnings`, not raised; an input that cannot be read at all (missing file, invalid gzip or UTF-8) still raises `InputFormatError`; invalid options raise `ConfigurationError`.

- **Example code:**

```python
import base64
import json

from dnakit.io_formats import validate_twobit, write_twobit

TWOBIT = base64.b64decode(
    write_twobit(">s1\nACGTNNacgt\n>s2\nGGCCAATT\n", verbose=False)["result"]["data_base64"]
)

envelope = validate_twobit(TWOBIT, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "valid": true,
    "errors": [],
    "warnings": [],
    "sequence_count": 2,
    "byte_order": "little",
    "file_size": 83
  }
}
```

## 122) validate_zarr

- **Purpose:** Validate a Zarr directory store and report its errors and warnings.

- **Method:** Opens the store, checks every array's shape, chunk shape and chunk counts, and decodes the first chunk of each initialised array so that corrupt codecs or truncated chunk files are detected (first chunks larger than `max_chunk_bytes` are not decoded).

- **API:** `dnakit.io_formats.validate_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: the path (`str` or `Path`) of a Zarr directory store (zarr v2 or v3; needs the optional `zarr>=3` and `numpy` packages); text content, bytes and gzip are not accepted. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_chunk_bytes`: int, default 67108864: first chunks larger than this are not decoded

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `valid`, `errors`, `warnings`, `array_count`, `group_count` and `path`. Content problems are reported in `errors`/`warnings`, not raised; an input that cannot be read at all (missing file, invalid gzip or UTF-8) still raises `InputFormatError`; invalid options raise `ConfigurationError`. If `zarr` (>= 3) or `numpy` cannot be imported the function returns `status: "not_implemented"` with `required_dependency: "pip install 'zarr>=3' numpy"` and a `message` (`zarr and numpy are required: <import error>`, or `zarr <version> is installed but this adapter needs zarr >= 3.`) instead of a `result`.

- **Example code:**

```python
import json

from dnakit.io_formats import validate_zarr, write_zarr

write_zarr(
    {"pos": [300, 100, 200, 400], "dosage": [[0, 1], [1, 1], [0, 0], [2, 1]]},
    options={"output": "demo.zarr", "overwrite": True},
    verbose=False,
)

envelope = validate_zarr("demo.zarr", verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "valid": true,
    "errors": [],
    "warnings": [],
    "array_count": 2,
    "group_count": 0,
    "path": "demo.zarr"
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

- **Purpose:** Serialise variants, or normalise existing `.map` content, into a PLINK `.map` file.

- **Method:** Every record is reduced to the four columns and written tab-separated, one variant per line, in the input order (cM with at most 15 significant digits). Existing `.map` text or a path is parsed first and then re-serialised, so whitespace-separated input comes out tab-separated. Nothing is written to disk unless `output` is given.

- **API:** `dnakit.io_formats.write_map(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a list of records, each a dict `{'chrom', 'id', 'cm', 'pos'}` or a list `[chrom, id, cm, bp]` / `[chrom, id, bp]` (`cm` then 0), or existing `.map` text/path (a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`) to normalise. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse text/file inputs larger than this
    - `output`: path to write (optional; the content is returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import write_map


envelope = write_map([["1", "rs1", 100], {"chrom": "2", "id": "rs2", "cm": 0.5, "pos": 10}], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "map",
    "record_count": 2,
    "bytes": 25,
    "sha256": "a0fc688e87258f9344dc06823654e51300a21c589a3421303fcada7d66454b49",
    "text": "1\trs1\t0\t100\n2\trs2\t0.5\t10\n"
  }
}
```

## 125) write_newick

- **Purpose:** Serialise trees, or normalise existing Newick content, into a Newick file.

- **Method:** Each tree is built from the given description (parsed Newick text, a flat node table, or a nested dict) and written as one normalised Newick line ending in `;`, in the input order. A flat table needs `id` and `parent` for each node (the root has `parent: null`); `name` and `length` are optional.

- **API:** `dnakit.io_formats.write_newick(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: Newick text/path (a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`), a `{'newick': ...}` dict, a flat `{'nodes': [{'id', 'parent', 'name', 'length'}, ...]}` dict, a nested `{'name', 'length', 'children'}` dict, or a list of these. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse text/file inputs larger than this
    - `output`: path to write (optional; the content is returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import write_newick


envelope = write_newick({"name": "root", "children": [{"name": "A", "length": 1}, {"name": "B", "length": 2}]}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "newick",
    "record_count": 1,
    "bytes": 15,
    "sha256": "d189445c127e5c47ea6623d600519888cb3f22fa38218bca1f679a1f5e099c38",
    "text": "(A:1,B:2)root;\n"
  }
}
```

## 126) write_nexus

- **Purpose:** Serialise sequences, or normalise existing NEXUS content, into a NEXUS file.

- **Method:** Writes a `#NEXUS` file with one DATA block (`DIMENSIONS NTAX/NCHAR`, `FORMAT DATATYPE=... MISSING=? GAP=-`, the matrix; names are quoted when needed), plus a TREES block when the input has trees. `layout='interleaved'` writes the matrix in blocks of `block_size` sites. FASTA, PHYLIP and NEXUS input is parsed first; plain mappings and record lists use `datatype` (default `dna`). Sequences of unequal length raise `InputFormatError`.

- **API:** `dnakit.io_formats.write_nexus(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: FASTA/PHYLIP/NEXUS text or path (a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`), a `{name: sequence}` mapping, or a list of `{'name', 'sequence'}` records. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse text/file inputs larger than this
    - `output`: path to write (optional; the content is returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file
    - `layout`: 'sequential' (default) or 'interleaved' for PHYLIP/NEXUS output
    - `block_size`: int, default 60: sites per block of an interleaved output
    - `strict_names`: bool, default False: write strict 10-character PHYLIP names (error when longer)
    - `line_width`: int, default 60: FASTA line width (0 = one line per sequence)
    - `datatype`: 'dna' (default), 'rna', 'protein' or 'any' for records without a declared type

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import write_nexus


envelope = write_nexus({"A": "ACGTACGT", "B": "ACGTACGA"}, options={"layout": "interleaved", "block_size": 4}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "nexus",
    "record_count": 2,
    "bytes": 178,
    "sha256": "90285f9c9dce3e7d66a73c016e18a2e474be23f9498a6210624311c0f37385f5",
    "text": "#NEXUS\n\nBEGIN DATA;\n    DIMENSIONS NTAX=2 NCHAR=8;\n    FORMAT DATATYPE=DNA MISSING=? GAP=- INTERLEAVE=YES;\n    MATRIX\n    A  ACGT\n    B  ACGT\n\n    A  ACGT\n    B  ACGA\n    ;\nEND;\n"
  }
}
```

## 127) write_paf

- **Purpose:** Serialise alignments, or normalise existing PAF content, into a PAF file.

- **Method:** Each record is written as one tab-separated line: the 12 mandatory fields followed by its tags as `TAG:TYPE:VALUE`. A tag given as a plain value in the input mapping (`{'tp': 'P'}`) is written with type `Z`. Existing PAF text or a path is parsed and re-serialised. The coordinates of each record are checked (`0 <= start <= end <= length` on query and target); an invalid record raises `InputFormatError`.

- **API:** `dnakit.io_formats.write_paf(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a list of dicts with the 12 mandatory PAF fields (`query_name`, `query_length`, `query_start`, `query_end`, `strand`, `target_name`, `target_length`, `target_start`, `target_end`, `matches`, `block_length`, `mapq`) and optional `tags`, or existing PAF text/path (a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`) to normalise. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse text/file inputs larger than this
    - `output`: path to write (optional; the content is returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import write_paf

RECORD = {
    "query_name": "q", "query_length": 10, "query_start": 0, "query_end": 10, "strand": "+",
    "target_name": "t", "target_length": 20, "target_start": 0, "target_end": 10,
    "matches": 9, "block_length": 10, "mapq": 5, "tags": {"tp": "P"},
}

envelope = write_paf([RECORD], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "paf",
    "record_count": 1,
    "bytes": 36,
    "sha256": "65ccfbff32154adca8949b2e2d4e75fb2fd8c528ca3f9274613f4f05b6df05a9",
    "text": "q\t10\t0\t10\t+\tt\t20\t0\t10\t9\t10\t5\ttp:Z:P\n"
  }
}
```

## 128) write_ped

- **Purpose:** Serialise individuals, or normalise existing `.ped` content, into a PLINK `.ped` file.

- **Method:** Each individual is written as `fid iid pat mat sex pheno` followed by its two alleles per variant, tab-separated, one individual per line, in the input order. Existing `.ped` text or a path is parsed first and re-serialised (whitespace-separated input comes out tab-separated). The genotype matrix is bounded by `max_cells`.

- **API:** `dnakit.io_formats.write_ped(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a list of individuals, each a dict `{'fid', 'iid', 'pat', 'mat', 'sex', 'pheno', 'genotypes': [[a1, a2], ...]}` (`fid` and `iid` are required; `pat`, `mat` and `sex` default to 0, `pheno` to -9 and `genotypes` to none), or existing `.ped` text/path (a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`) to normalise. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse text/file inputs larger than this
    - `output`: path to write (optional; the content is returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file
    - `max_cells`: int, default 20000000

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `variant_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import write_ped


envelope = write_ped([{"fid": "F1", "iid": "1", "sex": 1, "pheno": 2, "genotypes": [["A", "C"], ["G", "G"]]}], verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "ped",
    "record_count": 1,
    "variant_count": 2,
    "bytes": 21,
    "sha256": "df689e3f321aa8c782f823e51956fa08c4a41c841436038a2bde33b30e8f4944",
    "text": "F1\t1\t0\t0\t1\t2\tA\tC\tG\tG\n"
  }
}
```

## 129) write_phylip

- **Purpose:** Serialise sequences, or normalise existing PHYLIP content, into a PHYLIP alignment file.

- **Method:** Writes the header `ntax nchar` and one row per taxon: relaxed PHYLIP by default (name padded to 10 characters, then a space), or strict PHYLIP with `strict_names=True` (a name field of exactly 10 characters; a longer name raises `InputFormatError`). Names containing whitespace cannot be written. `layout='interleaved'` writes blocks of `block_size` sites. Sequences must have equal length (otherwise `InputFormatError`); the sequences are checked against the alphabet of `datatype` (`dna`, the default, and `rna`: IUPAC nucleotide codes, `-`, `?`, `.`, `*` and the gap/missing characters of the document, in any case; `protein`: the amino-acid letters plus B, Z, X, J, U, O, `*`, `-`, `?`, `.`; `any`: no check), and characters outside it raise `InputFormatError` (code `IO_FORMATS_BAD_ALPHABET`) that lists the offending characters and the sequences they occur in.

- **API:** `dnakit.io_formats.write_phylip(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: FASTA/PHYLIP/NEXUS text or path (a path (`str` or `Path`) or the content itself as `str`/`bytes` (a string that is an existing file is read as a path); gzip data is recognised by its magic bytes and unpacked, bounded by `max_bytes`), a `{name: sequence}` mapping, or a list of `{'name', 'sequence'}` records. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse text/file inputs larger than this
    - `output`: path to write (optional; the content is returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file
    - `datatype`: 'dna' (default), 'rna', 'protein' or 'any': alphabet the sequences are checked against (IUPAC codes, gap, ?, N for nucleotides; X for protein; 'any' skips the check)
    - `layout`: 'sequential' (default) or 'interleaved' for PHYLIP/NEXUS output
    - `block_size`: int, default 60: sites per block of an interleaved output
    - `strict_names`: bool, default False: write strict 10-character PHYLIP names (error when longer)
    - `line_width`: int, default 60: FASTA line width (0 = one line per sequence)
    - `input_layout`: 'auto' (default), 'sequential' or 'interleaved' when the input is PHYLIP text
    - `input_strict_names`: bool, default False when the input is PHYLIP text

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `bytes` and `sha256` of the serialised output, and the content as `text`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `text`. Malformed content or characters outside the alphabet raise `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import write_phylip


envelope = write_phylip(">a\nACGT\n>b\nACGA\n", options={"strict_names": True}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "phylip",
    "record_count": 2,
    "bytes": 35,
    "sha256": "162355e387d5b8a7ccc32c599fa049f4b8bdd217ad19b419ba040139e4875603",
    "text": " 2 4\na         ACGT\nb         ACGA\n"
  }
}
```

## 130) write_twobit

- **Purpose:** Serialise sequences, or re-encode an existing 2bit file, into a UCSC 2bit file.

- **Method:** Packs the sequences into a version-0 2bit file: A, C, G, T are stored 2 bits each, any other letter becomes an N block, and lower-case letters become mask blocks (so IUPAC codes such as R or Y are stored as N). `byte_order` selects little- (default) or big-endian. Sequence names must be unique and non-empty.

- **API:** `dnakit.io_formats.write_twobit(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: FASTA text/path, a `{name: sequence}` mapping, a list of `{'name', 'sequence'}` records, or an existing 2bit file (path or bytes) to re-encode. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `max_bytes`: int, default 268435456: refuse text/file inputs larger than this
    - `output`: path to write (optional; the content is returned inline otherwise)
    - `overwrite`: bool, default False: allow replacing an existing output file
    - `byte_order`: 'little' (default) or 'big': byte order of the written 2bit file

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count`, `sequence_count`, `bytes` and `sha256` of the serialised output, and the content as `data_base64`; with `options['output']` the file is written (an existing file needs `overwrite=True`) and `output_path` is returned instead of `data_base64`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import write_twobit


envelope = write_twobit(">s1\nACGTNNacgt\n>s2\nGGCC\n", verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "twobit",
    "record_count": 2,
    "sequence_count": 2,
    "bytes": 82,
    "sha256": "a75737f8defa5dad84385216e2f7b9091991fce11a1d9b3ebb45b4d15849423f",
    "data_base64": "QydBGgAAAAACAAAAAAAAAAJzMR4AAAACczJBAAAACgAAAAEAAAAEAAAAAgAAAAEAAAAGAAAABAAAAAAAAACcCcAEAAAAAAAAAAAAAAAAAAAA9Q=="
  }
}
```

## 131) write_zarr

- **Purpose:** Create a Zarr directory store from named arrays.

- **Method:** Creates a zarr v3 (default) or v2 (`zarr_format=2`) directory store at `output` with one array per name; `chunks` sets the chunk shape (default chosen by zarr) and `attrs` the attributes of the root group. An existing path needs `overwrite=True`; `output` is required. Arrays of other dtypes (for example strings) raise `InputFormatError`.

- **API:** `dnakit.io_formats.write_zarr(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a `{array name: array-like}` mapping, a list of `{'name', 'data'}` records, or one array-like (stored under `options['array']`, default `data`); supported dtypes are bool, int, uint and float. `options` is an optional mapping; unknown keys raise `ConfigurationError`. `verbose` (default `True`) prints progress messages; the examples pass `verbose=False`. Option keys (defaults as stated):

    - `output`: directory path of the new zarr store (required)
    - `overwrite`: bool, default False: replace an existing zarr store at output
    - `array`: name of the array when the input is a single array-like (default 'data')
    - `chunks`: int or list of ints: chunk shape (default: chosen by zarr)
    - `attrs`: mapping stored as attributes of the root group
    - `zarr_format`: 2 or 3 (default 3)

- **Output:** The call returns the DNAKit envelope with `status` (`"ok"`), `module` (`"io_formats"`), `operation`, `provenance` (`implementation: "native"`, `backend: "local"`), `input_summary`, `warnings` and `result` (the example prints `status`, `warnings` and `result`). `result`: `format`, `record_count` (number of arrays written), `zarr_format` and `output_path` (the store is written to `options['output']`; no content is returned). If `zarr` (>= 3) or `numpy` cannot be imported the function returns `status: "not_implemented"` with `required_dependency: "pip install 'zarr>=3' numpy"` and a `message` (`zarr and numpy are required: <import error>`, or `zarr <version> is installed but this adapter needs zarr >= 3.`) instead of a `result`. Malformed content raises `InputFormatError`; invalid or unknown options, and inputs above the size limits, raise `ConfigurationError`.

- **Example code:**

```python
import json

from dnakit.io_formats import write_zarr


envelope = write_zarr({"pos": [100, 200, 300]}, options={"output": "written.zarr", "overwrite": True, "chunks": 2, "attrs": {"source": "demo"}}, verbose=False)
print(json.dumps({key: envelope[key] for key in ("status", "warnings", "result")}, indent=2))
```

- **Example output:**

```json
{
  "status": "ok",
  "warnings": [],
  "result": {
    "format": "zarr",
    "record_count": 1,
    "zarr_format": 3,
    "output_path": "written.zarr"
  }
}
```

