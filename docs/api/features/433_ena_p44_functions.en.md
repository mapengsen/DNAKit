# P44 Public database operations

Execute real queries and downloads against fixed official endpoints. Ordinary tests and the replay example below use no network; they do not claim per-record live validation.

Downloads inherit byte limits, atomic writes, checksums and progress. Results are JSON serialisable and retain responses, retrieval dates and page information. Follow upstream page tokens explicitly; list release files before choosing an archived dataset.

## `ena_count_keyword_results()`

- **Purpose:** Count ENA keyword-search hits for a record format and domain.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ena_count_keyword_results(query: 'str', domain: 'str', result: 'str', *, format: 'str' = 'embl', verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ena_download_cram_reference()`

- **Purpose:** Download exact reference sequence bytes by CRAM registry MD5 or SHA1 checksum.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ena_download_cram_reference(checksum: 'str', output_path: 'str', *, algorithm: 'str' = 'md5', overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ena_download_keyword_results()`

- **Purpose:** Download keyword-search EMBL/FASTA/XML records, including WGS/TSA/TLS master sets.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ena_download_keyword_results(query: 'str', domain: 'str', result: 'str', output_path: 'str', *, format: 'str' = 'embl', limit: 'int' = 100000, offset: 'int' = 0, overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ena_download_linked_records()`

- **Purpose:** Download EMBL/FASTA/XML records linked to a study, sample or taxon accession.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ena_download_linked_records(link_type: 'str', accession: 'str', result: 'str', output_path: 'str', *, format: 'str' = 'embl', subtree: 'bool' = False, overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ena_download_record()`

- **Purpose:** Download EMBL contig-set records or XML assembly records by accession.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ena_download_record(accession: 'str', output_path: 'str', *, format: 'str' = 'embl', full_contig_set: 'bool' = False, overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ena_download_search()`

- **Purpose:** Download advanced-search EMBL/FASTA sequences or XML run/experiment/analysis/assembly
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ena_download_search(query: 'str', result: 'str', output_path: 'str', *, format: 'str' = 'embl', limit: 'int' = 100000, overwrite: 'bool' = False, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ena_linked_accessions()`

- **Purpose:** List portal DNA record accessions linked to a study, sample or taxon.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ena_linked_accessions(link_type: 'str', accession: 'str', result: 'str', *, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ena_search_cross_references()`

- **Purpose:** Search ENA cross-references to an external database and return original TSV.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ena_search_cross_references(source: 'str', *, source_accession: 'str | None' = None, target: 'str | None' = None, target_accession: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## Offline input validation example

**Example code:**

```python
from dnakit.download.ena_p44_functions import ena_download_record
from dnakit.exceptions import ConfigurationError
try:
    ena_download_record("X", "unused", format="bad", verbose=False)
except ConfigurationError:
    print("invalid input rejected before network access")
```

**Example output:**

```text
invalid input rejected before network access
```

## Upstream references

- [官方文档 1](https://rest.ensembl.org/)
- [官方文档 2](https://github.com/ncbi/datasets/blob/master/datasets.openapi.yaml)
- [官方文档 3](https://www.ebi.ac.uk/ena/browser/api/swagger-ui/index.html)
