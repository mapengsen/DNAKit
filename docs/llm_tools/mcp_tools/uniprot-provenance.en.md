# UniProt MCP Server (smaniches) (`uniprot-provenance`)

DNAKit's MCP tool index keeps **7** DNA/gene tools of this server (41 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [smaniches/uniprot-mcp @ `a43171f`](https://github.com/smaniches/uniprot-mcp/tree/a43171f41a8c9119883834ceb02b2eb1fa64b415) |
| License | Apache-2.0 |
| Command | `uniprot-mcp-server` |
| Transports | stdio |
| Installation | pip install uniprot-mcp-server (PyPI 1.3.9) into a Python &gt;=3.11 environment; optionally set UNIPROT_PIN_RELEASE (for example 2026_01) to pin the UniProt release and UNIPROT_MCP_CACHE_DIR to enable the offline replay cache directory. No API key is needed. |
| Effects | Read-only HTTPS queries to rest.uniprot.org and www.uniprot.org, NCBI E-utilities (ClinVar) and the AlphaFold EBI service; ID mapping submits identifiers to UniProt as a job, and uniprot_replay_from_cache reads local files from UNIPROT_MCP_CACHE_DIR; live responses are not written to that cache in the pinned code. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`uniprot_get_disease_associations`](#uniprot_get_disease_associations) | `association` | Returns UniProt-curated disease associations for an entry with names, acronyms, disease IDs and OMIM cross-references. |
| [`uniprot_get_go_terms`](#uniprot_get_go_terms) | `genome_annotation` | Returns Gene Ontology annotations of a UniProt entry grouped by aspect (function, process, component). |
| [`uniprot_get_variants`](#uniprot_get_variants) | `variant` | Lists literature-described natural variants of a UniProt entry, including disease-associated mutations. |
| [`uniprot_id_mapping`](#uniprot_id_mapping) | `genome_annotation` | Maps up to 100 identifiers between UniProt and external databases such as Ensembl, Entrez GeneID or PDB. |
| [`uniprot_lookup_variant`](#uniprot_lookup_variant) | `variant` | Looks up a protein change (e.g. R175H) in a UniProt entry's natural-variant annotations with curated description. |
| [`uniprot_resolve_clinvar`](#uniprot_resolve_clinvar) | `variant` | Looks up ClinVar records for the gene of a UniProt entry, optionally filtered by protein change, with significance. |
| [`uniprot_resolve_orthology`](#uniprot_resolve_orthology) | `genome_annotation` | Groups a UniProt entry's orthology cross-references by database (KEGG, OMA, OrthoDB, eggNOG, PANTHER, GeneTree). |

## `uniprot_get_disease_associations` {#uniprot_get_disease_associations}

Returns UniProt-curated disease associations for an entry with names, acronyms, disease IDs and OMIM cross-references.

- **Category:** `association` (GWAS, QTL, LD and gene/variant-disease or trait associations)
- **Required arguments:** `accession`
- **Keywords:** disease association, omim, mim, curated, gene-disease
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | UniProt accession, e.g. 'P04637' (human TP53) or 'P38398' (human BRCA1). Both reviewed (Swiss-Prot) and unreviewed (TrEMBL) accessions are accepted; case-sensitive. |
| `response_format` | string | no | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. Default `markdown` |

??? note "Upstream description"

    Structured disease associations for a UniProt entry. Returns the
    diseases recorded in DISEASE-type comments with name, acronym,
    UniProt disease ID, OMIM cross-reference, description, and the
    annotation note. Critical for clinical interpretation — distinguishes
    a UniProt-curated disease association (literature-anchored) from a
    raw cross-reference. Empty result does not imply disease-irrelevant;
    see Open Targets / OMIM / DisGeNET for population-level evidence.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'uniprot-provenance',
    'uniprot_get_disease_associations',
    {'accession': '...'},
)
```

## `uniprot_get_go_terms` {#uniprot_get_go_terms}

Returns Gene Ontology annotations of a UniProt entry grouped by aspect (function, process, component).

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** gene ontology, go, molecular function, biological process, cellular component
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | UniProt accession, e.g. 'P04637' (human TP53) or 'P38398' (human BRCA1). Both reviewed (Swiss-Prot) and unreviewed (TrEMBL) accessions are accepted; case-sensitive. |
| `aspect` | string | no | Optional Gene Ontology aspect filter: 'F' (molecular function), 'P' (biological process), 'C' (cellular component), or empty for all three. Default `` |
| `response_format` | string | no | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. Default `markdown` |

??? note "Upstream description"

    Get GO annotations grouped by aspect.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'uniprot-provenance',
    'uniprot_get_go_terms',
    {'accession': '...'},
)
```

## `uniprot_get_variants` {#uniprot_get_variants}

Lists literature-described natural variants of a UniProt entry, including disease-associated mutations.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `accession`
- **Keywords:** natural variant, missense, disease mutation, amino acid change
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | UniProt accession, e.g. 'P04637' (human TP53) or 'P38398' (human BRCA1). Both reviewed (Swiss-Prot) and unreviewed (TrEMBL) accessions are accepted; case-sensitive. |
| `response_format` | string | no | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. Default `markdown` |

??? note "Upstream description"

    List every literature-described natural variant UniProt has curated
    for an entry, including disease-associated mutations. Use this to see
    the full variant catalogue for a protein; to check one specific
    HGVS-shorthand change (e.g. 'R175H') use ``uniprot_lookup_variant``
    instead, which does the position/residue matching for you. UniProt's
    natural-variant annotations only cover literature-described variants —
    for population-scale clinical significance data use
    ``uniprot_resolve_clinvar``.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'uniprot-provenance',
    'uniprot_get_variants',
    {'accession': '...'},
)
```

## `uniprot_id_mapping` {#uniprot_id_mapping}

Maps up to 100 identifiers between UniProt and external databases such as Ensembl, Entrez GeneID or PDB.

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `ids`, `from_db`, `to_db`
- **Keywords:** id mapping, ensembl, entrez, geneid, cross-reference
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ids` | string | yes | Comma-separated identifiers to map, up to 100 per call. |
| `from_db` | string | yes | Source database code, e.g. 'UniProtKB_AC-ID', 'PDB', 'Ensembl', 'GeneID' (Entrez), or 'Gene_Name'. |
| `to_db` | string | yes | Target database code, same code set as ``from_db``. |
| `response_format` | string | no | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. Default `markdown` |

??? note "Upstream description"

    Map identifiers between UniProt and external databases (or between
    two external databases) via UniProt's ID mapping service. Submits an
    async job and polls it to completion server-side, so the call may take
    a few seconds for large batches.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'uniprot-provenance',
    'uniprot_id_mapping',
    {'ids': '...', 'from_db': '...', 'to_db': '...'},
)
```

## `uniprot_lookup_variant` {#uniprot_lookup_variant}

Looks up a protein change (e.g. R175H) in a UniProt entry's natural-variant annotations with curated description.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `accession`, `change`
- **Keywords:** amino acid change, r175h, natural variant, disease association, hgvs shorthand
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | UniProt accession, e.g. 'P04637' (human TP53) or 'P38398' (human BRCA1). Both reviewed (Swiss-Prot) and unreviewed (TrEMBL) accessions are accepted; case-sensitive. |
| `change` | string | yes | HGVS-shorthand amino-acid change, e.g. 'R175H', 'V600E', 'R248*' (stop). Format: &lt;original residue&gt;&lt;1-indexed position&gt;&lt;alt residue or '*'&gt;. |
| `response_format` | string | no | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. Default `markdown` |

??? note "Upstream description"

    Look up an HGVS-shorthand amino-acid change (e.g. ``R175H``,
    ``V600E``, ``R248*``) in the UniProt entry's natural-variant
    annotations. Returns the matching variant feature(s) including the
    UniProt-curated description (often a disease association). A null
    result here does NOT mean a variant is benign — UniProt only
    annotates literature-described variants; ClinVar / dbSNP carry
    population-level data.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'uniprot-provenance',
    'uniprot_lookup_variant',
    {'accession': '...', 'change': '...'},
)
```

## `uniprot_resolve_clinvar` {#uniprot_resolve_clinvar}

Looks up ClinVar records for the gene of a UniProt entry, optionally filtered by protein change, with significance.

- **Category:** `variant` (Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores)
- **Required arguments:** `accession`
- **Keywords:** clinvar, clinical significance, review status, protein change, gene symbol
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | UniProt accession, e.g. 'P04637' (human TP53) or 'P38398' (human BRCA1). Both reviewed (Swiss-Prot) and unreviewed (TrEMBL) accessions are accepted; case-sensitive. |
| `change` | string | no | Optional HGVS-shorthand protein change to filter to, e.g. 'R175H'. Leave empty to return all ClinVar records for the gene. Default `` |
| `size` | integer | no | Maximum number of ClinVar records to return; capped at 50. Default `10` |
| `response_format` | string | no | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. Default `markdown` |

??? note "Upstream description"

    Look up ClinVar records for the gene encoded by a UniProt entry.
    First fetches the entry to extract the canonical gene symbol, then
    queries NCBI eutils ClinVar by gene (and optional protein-change
    filter, e.g. ``R175H``). Returns clinical-significance classification,
    review status, condition list (trait_set), molecular consequence,
    and the protein-change list per record.

    Critical for clinical workflows — UniProt's natural-variant
    annotations stop at literature-described variants. ClinVar carries
    every variant submitted by clinical labs, with curated significance
    classifications. Combine ``uniprot_lookup_variant`` (UniProt side)
    with ``uniprot_resolve_clinvar`` (population side) for a full
    variant-effect picture.

    Calls https://eutils.ncbi.nlm.nih.gov — declared in PRIVACY.md.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'uniprot-provenance',
    'uniprot_resolve_clinvar',
    {'accession': '...'},
)
```

## `uniprot_resolve_orthology` {#uniprot_resolve_orthology}

Groups a UniProt entry's orthology cross-references by database (KEGG, OMA, OrthoDB, eggNOG, PANTHER, GeneTree).

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `accession`
- **Keywords:** orthologs, orthology, eggnog, oma, panther, cross-species
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | UniProt accession, e.g. 'P04637' (human TP53) or 'P38398' (human BRCA1). Both reviewed (Swiss-Prot) and unreviewed (TrEMBL) accessions are accepted; case-sensitive. |
| `response_format` | string | no | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. Default `markdown` |

??? note "Upstream description"

    Group every orthology cross-reference in a UniProt entry by source
    database (KEGG / OMA / OrthoDB / eggNOG / HOGENOM / PhylomeDB /
    InParanoid / TreeFam / GeneTree / PAN-GO / PANTHER / OrthoInspector).
    Different databases use different inference methods; surfacing them
    side-by-side lets the agent reason about consensus when comparing
    orthologs across species. Pure-Python — no extra HTTP call beyond
    the entry fetch.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'uniprot-provenance',
    'uniprot_resolve_orthology',
    {'accession': '...'},
)
```
