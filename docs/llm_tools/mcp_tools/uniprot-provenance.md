# UniProt MCP Server (smaniches) (`uniprot-provenance`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **7** 个 DNA/基因相关工具（服务器共列出 41 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [smaniches/uniprot-mcp @ `a43171f`](https://github.com/smaniches/uniprot-mcp/tree/a43171f41a8c9119883834ceb02b2eb1fa64b415) |
| 许可证 | Apache-2.0 |
| 启动命令 | `uniprot-mcp-server` |
| 传输方式 | stdio |
| 安装 | pip install uniprot-mcp-server (PyPI 1.3.9) into a Python &gt;=3.11 environment; optionally set UNIPROT_PIN_RELEASE (for example 2026_01) to pin the UniProt release and UNIPROT_MCP_CACHE_DIR to enable the offline replay cache directory. No API key is needed. |
| 副作用 | Read-only HTTPS queries to rest.uniprot.org and www.uniprot.org, NCBI E-utilities (ClinVar) and the AlphaFold EBI service; ID mapping submits identifiers to UniProt as a job, and uniprot_replay_from_cache reads local files from UNIPROT_MCP_CACHE_DIR; live responses are not written to that cache in the pinned code. |

| 工具 | 类别 | 简介 |
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

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`accession`
- **关键词：**disease association, omim, mim, curated, gene-disease
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | UniProt accession, e.g. 'P04637' (human TP53) or 'P38398' (human BRCA1). Both reviewed (Swiss-Prot) and unreviewed (TrEMBL) accessions are accepted; case-sensitive. |
| `response_format` | string | 否 | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. 默认 `markdown` |

??? note "上游说明（英文原文）"

    Structured disease associations for a UniProt entry. Returns the
    diseases recorded in DISEASE-type comments with name, acronym,
    UniProt disease ID, OMIM cross-reference, description, and the
    annotation note. Critical for clinical interpretation — distinguishes
    a UniProt-curated disease association (literature-anchored) from a
    raw cross-reference. Empty result does not imply disease-irrelevant;
    see Open Targets / OMIM / DisGeNET for population-level evidence.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`accession`
- **关键词：**gene ontology, go, molecular function, biological process, cellular component
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | UniProt accession, e.g. 'P04637' (human TP53) or 'P38398' (human BRCA1). Both reviewed (Swiss-Prot) and unreviewed (TrEMBL) accessions are accepted; case-sensitive. |
| `aspect` | string | 否 | Optional Gene Ontology aspect filter: 'F' (molecular function), 'P' (biological process), 'C' (cellular component), or empty for all three. 默认 `` |
| `response_format` | string | 否 | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. 默认 `markdown` |

??? note "上游说明（英文原文）"

    Get GO annotations grouped by aspect.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`accession`
- **关键词：**natural variant, missense, disease mutation, amino acid change
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | UniProt accession, e.g. 'P04637' (human TP53) or 'P38398' (human BRCA1). Both reviewed (Swiss-Prot) and unreviewed (TrEMBL) accessions are accepted; case-sensitive. |
| `response_format` | string | 否 | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. 默认 `markdown` |

??? note "上游说明（英文原文）"

    List every literature-described natural variant UniProt has curated
    for an entry, including disease-associated mutations. Use this to see
    the full variant catalogue for a protein; to check one specific
    HGVS-shorthand change (e.g. 'R175H') use ``uniprot_lookup_variant``
    instead, which does the position/residue matching for you. UniProt's
    natural-variant annotations only cover literature-described variants —
    for population-scale clinical significance data use
    ``uniprot_resolve_clinvar``.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`ids`, `from_db`, `to_db`
- **关键词：**id mapping, ensembl, entrez, geneid, cross-reference
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ids` | string | 是 | Comma-separated identifiers to map, up to 100 per call. |
| `from_db` | string | 是 | Source database code, e.g. 'UniProtKB_AC-ID', 'PDB', 'Ensembl', 'GeneID' (Entrez), or 'Gene_Name'. |
| `to_db` | string | 是 | Target database code, same code set as ``from_db``. |
| `response_format` | string | 否 | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. 默认 `markdown` |

??? note "上游说明（英文原文）"

    Map identifiers between UniProt and external databases (or between
    two external databases) via UniProt's ID mapping service. Submits an
    async job and polls it to completion server-side, so the call may take
    a few seconds for large batches.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`accession`, `change`
- **关键词：**amino acid change, r175h, natural variant, disease association, hgvs shorthand
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | UniProt accession, e.g. 'P04637' (human TP53) or 'P38398' (human BRCA1). Both reviewed (Swiss-Prot) and unreviewed (TrEMBL) accessions are accepted; case-sensitive. |
| `change` | string | 是 | HGVS-shorthand amino-acid change, e.g. 'R175H', 'V600E', 'R248*' (stop). Format: &lt;original residue&gt;&lt;1-indexed position&gt;&lt;alt residue or '*'&gt;. |
| `response_format` | string | 否 | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. 默认 `markdown` |

??? note "上游说明（英文原文）"

    Look up an HGVS-shorthand amino-acid change (e.g. ``R175H``,
    ``V600E``, ``R248*``) in the UniProt entry's natural-variant
    annotations. Returns the matching variant feature(s) including the
    UniProt-curated description (often a disease association). A null
    result here does NOT mean a variant is benign — UniProt only
    annotates literature-described variants; ClinVar / dbSNP carry
    population-level data.

安装服务器后通过 DNAKit 调用：

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

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`accession`
- **关键词：**clinvar, clinical significance, review status, protein change, gene symbol
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | UniProt accession, e.g. 'P04637' (human TP53) or 'P38398' (human BRCA1). Both reviewed (Swiss-Prot) and unreviewed (TrEMBL) accessions are accepted; case-sensitive. |
| `change` | string | 否 | Optional HGVS-shorthand protein change to filter to, e.g. 'R175H'. Leave empty to return all ClinVar records for the gene. 默认 `` |
| `size` | integer | 否 | Maximum number of ClinVar records to return; capped at 50. 默认 `10` |
| `response_format` | string | 否 | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. 默认 `markdown` |

??? note "上游说明（英文原文）"

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

安装服务器后通过 DNAKit 调用：

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

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`accession`
- **关键词：**orthologs, orthology, eggnog, oma, panther, cross-species
- **只读：**是（服务器声明）

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | UniProt accession, e.g. 'P04637' (human TP53) or 'P38398' (human BRCA1). Both reviewed (Swiss-Prot) and unreviewed (TrEMBL) accessions are accepted; case-sensitive. |
| `response_format` | string | 否 | 'markdown' (default) for a human-readable report with a provenance footer, or 'json' for a machine-parseable structured payload with the same data. Any other value is rejected. 默认 `markdown` |

??? note "上游说明（英文原文）"

    Group every orthology cross-reference in a UniProt entry by source
    database (KEGG / OMA / OrthoDB / eggNOG / HOGENOM / PhylomeDB /
    InParanoid / TreeFam / GeneTree / PAN-GO / PANTHER / OrthoInspector).
    Different databases use different inference methods; surfacing them
    side-by-side lets the agent reason about consensus when comparing
    orthologs across species. Pure-Python — no extra HTTP call beyond
    the entry fetch.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'uniprot-provenance',
    'uniprot_resolve_orthology',
    {'accession': '...'},
)
```
