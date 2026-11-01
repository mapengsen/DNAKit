# P43 Public database operations

Execute real queries and downloads against fixed official endpoints. Ordinary tests and the replay example below use no network; they do not claim per-record live validation.

Downloads inherit byte limits, atomic writes, checksums and progress. Results are JSON serialisable and retain responses, retrieval dates and page information. Follow upstream page tokens explicitly; list release files before choosing an archived dataset.

## `ncbi_fetch_dna_records()`

- **Purpose:** Fetch nucleotide FASTA/gbc/ft/native/CDS, Gene XML/table/ASN.1, SNP, ClinVar or SRA records.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ncbi_fetch_dna_records(database: 'str', ids: 'Sequence[str]', *, rettype: 'str' = 'fasta', retmode: 'str' = 'text', email: 'str | None' = None, api_key: 'str | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ncbi_link_dna_records()`

- **Purpose:** ELink one of the 28 curated nuccore/assembly/gene/variant/SRA/GEO/OMIM relationships.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ncbi_link_dna_records(linkname: 'str', ids: 'Sequence[str]', *, email: 'str | None' = None, api_key: 'str | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ncbi_search_dna_database()`

- **Purpose:** ESearch assembly, annotinfo, snp, dbvar, grasp, geoprofiles, gtr or seqannot; return IDs
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ncbi_search_dna_database(database: 'str', term: 'str', *, retmax: 'int' = 20, retstart: 'int' = 0, email: 'str | None' = None, api_key: 'str | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## `ncbi_summarize_dna_records()`

- **Purpose:** Return unmodified ESummary JSON records from snp, sra or omim.
- **Method:** Call its upstream API; download functions write actual bytes and raise on failure.
- **API:** `ncbi_summarize_dna_records(database: 'str', ids: 'Sequence[str]', *, email: 'str | None' = None, api_key: 'str | None' = None, timeout: 'float' = 60, max_bytes: 'int' = 20000000, verbose: 'bool' = True) -> 'dict[str, Any]'`
- **Input:** Required identifiers, collection selectors and optional parameters are shown in the signature; names follow upstream documentation.
- **Output:** `status`, `data` or downloaded-file metadata, plus provenance.

## Offline input validation example

**Example code:**

```python
from dnakit.databases.ncbi_eutils_p43_functions import ncbi_fetch_dna_records
from dnakit.exceptions import ConfigurationError
try:
    ncbi_fetch_dna_records("protein", ["1"], verbose=False)
except ConfigurationError:
    print("invalid input rejected before network access")
```

**Example output:**

```text
invalid input rejected before network access
```

## Upstream references

- [官方文档 1](https://www.ncbi.nlm.nih.gov/books/NBK25499/)
