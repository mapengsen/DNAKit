---
name: sequence-retrieval
description: "Retrieve DNA, RNA and protein sequences with DNAKit from NCBI and ENA. Use for accession/version lookup, organism and gene disambiguation, transcript isoforms, reference-versus-submission comparisons, and FASTA/GenBank/EMBL downloads with source records."
disable-model-invocation: true
---

## DNAKit execution

Use `dnakit skills doctor sequence-retrieval` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# Biological Sequence Retrieval

Use DNAKit to identify the intended sequence, retrieve it, and produce a sequence
profile. This adaptation retains the upstream workflow; see [UPSTREAM.md](UPSTREAM.md)
for authorship, license and changes. ToolUniverse is no longer a runtime dependency.

## Dependencies and entry points

The base DNAKit package performs NCBI search, UID resolution, FASTA and GenBank
downloads, and ENA metadata/FASTA/EMBL retrieval. Public services require network
access. Only format conversion requires `pip install 'dnakit[bio-formats]'`:
NCBI GenBank → EMBL, or ENA EMBL → GenBank. Biopython handles those conversions.

Python uses `import dnakit`. The equivalent tools exposed by DNAKit MCP are
`dnakit.search.entrez`, `dnakit.search.accessions`, `dnakit.search.accession`,
`dnakit.search.ena_search`, `dnakit.download.ncbi_sequence` and
`dnakit.download.ena_sequence`. Describe a tool before constructing an MCP call;
use its actual input schema and the client's normal file-write authorization.

## 1. Clarify and disambiguate

Identify organism, gene or locus, molecule (genomic DNA / transcript / protein),
strain when relevant, and output format. Ask when the answer changes which
sequence is appropriate; a specific accession normally needs no clarification.
Search scientific names and established gene symbols. Report in the user's language.

Verify accessions and versions from provider responses. Do not invent an accession,
choose a canonical isoform solely from search ranking, or silently replace a
requested version. For human transcript work, inspect MANE Select evidence when
available, and retain the user's intended transcript or isoform.

| Identifier | Provider selection |
| --- | --- |
| RefSeq, including NC_, NG_, NM_, NR_, NP_, NT_, NW_, NZ_, XM_, XR_, XP_, WP_, YP_ | NCBI; do not pass the RefSeq identifier to ENA |
| INSDC nucleotide accessions such as U00096.3 or CP… | NCBI or ENA, subject to record availability |
| GCA_/GCF_ assembly accessions | Use DNAKit genome-package download, not nucleotide EFetch |
| Protein accession | NCBI with `database="protein"` |

RefSeq status, predicted status, annotation evidence and relevance are separate
properties. XM_/XR_/XP_ indicate predicted records. A prefix alone does not
establish experimental validation or a universal sequence-quality ranking.

## 2. Search and retrieve

```python
import dnakit

result = dnakit.search.ncbi.entrez(
    "nuccore", organism="Escherichia coli", strain="K-12",
    seq_type="complete_genome", retmax=3,
)
accessions = dnakit.search.ncbi.accessions(result.metadata["uids"])
print(result.total_count, accessions)
for summary in result.records:
    print(summary)  # Inspect title, organism, accession version and length.
```

`entrez` retains positional `database, term` support. Structured filters are ANDed
with an explicit term. `organism`, `gene`, and `strain` are literal values;
`keywords` and `term` accept Entrez expressions. Strain uses All Fields, so inspect
the actual strain in the resulting annotation. Available `seq_type` values are
`complete_genome`, `genomic`, `mrna`, `ncrna`, `refseq`, and `protein`; protein
requires `database="protein"`. `source="refseq"` can be combined with a molecule
type; `source="genbank"` excludes RefSeq without asserting curation quality.

`retmax` bounds this page and `retstart` selects later pages. `total_count` is the
total search count; `records` holds summaries and `metadata["uids"]` holds UIDs.
`accessions` returns versioned strings in UID input order, preserves duplicates,
and returns an empty tuple for empty input. Do not index empty results.

After selecting and verifying the intended accession:

```python
download = dnakit.download.ncbi.sequence(
    "NC_000913.3", "ecoli_k12.fasta", format="fasta",
)
print(download.files[0].path, download.manifest_path)
```

This downloads a fixed reference, not an automatically selected search hit.
Use the accession chosen from the summaries when the task concerns search results.
Downloads return `DatasetDownloadResult` and save a sibling `.manifest.json`
containing the requested accession, URL, format/conversion, byte size and hashes.
An unversioned request may resolve to a versioned record: report the identifier
in the actual file. Existing files are protected unless `DownloadConfig(overwrite=True)`
is explicitly appropriate. `len(fasta_text)` includes headers/newlines; count
sequence symbols or use record metadata for biological length.

Fetch `format="genbank"` for feature annotations (GenPept for protein). Nucleotide
`format="embl"` converts the GenBank response; it does not imply an ENA accession.
Optional `seq_start`, `seq_stop` use **1-based inclusive** coordinates and must be
provided together; `strand=2` requests the reverse complement. Do not mix these
with DNAKit's other APIs that use 0-based half-open coordinates.

For INSDC nucleotide records, use `dnakit.search.ena.ena_search` with
`result="sequence"` for metadata, and `dnakit.download.ena.sequence` for files.
Use [examples.md](examples.md) for complete independent examples, including ENA,
RNA/protein retrieval and multiple formats.

## 3. Handle errors and reconcile sources

- No hits: inspect spelling, synonyms and filters before broadening the query.
  Distinguish no results from network, service or invalid-query failures.
- NCBI unavailable: ENA is an alternative for the same INSDC accession. A RefSeq
  → INSDC change requires a verified cross-reference and an explicit report of it.
- ENA 404: check availability and identifier type; 404 alone does not prove RefSeq.
- Conflicting records: compare version, strain, coverage and annotation before
  interpreting discrepancies. Do not assume that every difference is an error correction.
- Responses are bounded (20 MB by default, at most 100 MB via `SearchConfig`).
  For larger genomes use assembly-package downloads or provide the exact provider
  link and explain what was not retrieved. Avoid dumping full genomes into context.
- `SearchConfig` accepts optional NCBI `api_key` and `email`; pass the same config
  to search/resolution and as `search_config` to sequence downloads. Never put
  credentials into example files or reports. Provider failures are retried under
  DNAKit's bounded network retry policy; bulk jobs must respect provider rate limits.

## 4. Report a sequence profile

Use [CHECKLIST.md](CHECKLIST.md) to verify completeness. Include query, database,
total hits, chosen accession/version, organism/strain, molecule, biological length,
topology, reference/prediction status and evidence. Show a short FASTA preview and
the saved file/manifest paths. From annotated files report CDS, gene, tRNA, rRNA
and regulatory-feature counts where present; mark unavailable fields as unknown.

List relevant alternative isoforms/sequences and explain their relevance. Include
verified RefSeq/INSDC, BioProject and BioSample cross-references; do not infer links
merely from accession similarity. Give format/download options, retrieval date,
source links and citations to DNAKit, the upstream Skill and the actual databases.
Summarize decision-relevant limitations without narrating every request.
