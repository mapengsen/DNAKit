# P45 Public database operations

Execute real queries and downloads against fixed official endpoints. Ordinary tests and the replay example below use no network; they do not claim per-record live validation.

Downloads inherit byte limits, atomic writes, checksums and progress. Results are JSON serialisable and retain responses, retrieval dates and page information. Follow upstream page tokens explicitly; list release files before choosing an archived dataset.

## `encode_body_system_summary()`

- **Purpose:** Aggregate ENCODE Experiment counts by body system and organ.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `encode_body_system_summary(*, filters: 'Mapping[str, Any] | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `encode_download_urls()`

- **Purpose:** Retrieve ENCODE batch-download files.txt URL list for a search or cart.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `encode_download_urls(*, filters: 'Mapping[str, Any] | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `encode_experiment_matrix()`

- **Purpose:** Aggregate ENCODE Experiment counts in the assay x biosample matrix.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `encode_experiment_matrix(*, filters: 'Mapping[str, Any] | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `encode_object_metadata()`

- **Purpose:** Retrieve ENCODE object metadata by accession or @id, using frame=object.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `encode_object_metadata(accession: 'str', *, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `encode_tabular_report()`

- **Purpose:** Return ENCODE report.tsv for selected fields and a bounded or explicit all-results search.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `encode_tabular_report(fields: 'Sequence[str]', *, filters: 'Mapping[str, Any] | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `rcsb_nucleotide_motif_search()`

- **Purpose:** Search PDB DNA/RNA chain motifs in simple, prosite or regex notation.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `rcsb_nucleotide_motif_search(pattern: 'str', *, sequence_type: 'str' = 'dna', pattern_type: 'str' = 'simple', start: 'int' = 0, rows: 'int' = 10, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `rcsb_nucleotide_sequence_search()`

- **Purpose:** Search PDB DNA or RNA chains by sequence similarity, returning polymer-entity hits.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `rcsb_nucleotide_sequence_search(sequence: 'str', *, sequence_type: 'str' = 'dna', identity_cutoff: 'float' = 0.0, evalue_cutoff: 'float' = 0.1, start: 'int' = 0, rows: 'int' = 10, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ucsc_search_assembly()`

- **Purpose:** Search UCSC track descriptions, names and genomic feature data in one assembly.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ucsc_search_assembly(genome: 'str', term: 'str', *, categories: 'str | None' = None, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ucsc_track_schema()`

- **Purpose:** Return original UCSC annotation-track column schema.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ucsc_track_schema(genome: 'str', track: 'str', *, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## Offline input validation example

**Example code:**

```python
from dnakit.databases.public_apis_p45_functions import ucsc_track_schema
from dnakit.exceptions import ConfigurationError
try:
    ucsc_track_schema("", "knownGene", verbose=False)
except ConfigurationError:
    print("invalid input rejected before network access")
```

**Example output:**

```text
invalid input rejected before network access
```

## Upstream references

- [官方文档 1](https://genome.ucsc.edu/goldenPath/help/api.html)
- [官方文档 2](https://www.encodeproject.org/help/rest-api/)
- [官方文档 3](https://search.rcsb.org/)
