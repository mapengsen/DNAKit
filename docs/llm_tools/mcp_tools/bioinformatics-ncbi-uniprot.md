# Bioinformatics MCP Server (NCBI and UniProt) (`bioinformatics-ncbi-uniprot`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **4** 个 DNA/基因相关工具（服务器共列出 10 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [AiAgentKarl/bioinformatics-mcp-server @ `b073742`](https://github.com/AiAgentKarl/bioinformatics-mcp-server/tree/b073742a7059f5cad816d0cf41f9fcdf49586933) |
| 许可证 | MIT |
| 启动命令 | `bioinformatics-server` |
| 传输方式 | stdio |
| 安装 | pip install bioinformatics-mcp-server (PyPI 0.1.0) into a Python &gt;=3.11 environment; no environment variables or API keys are needed. The wheel installs a top-level package named src, so use a dedicated environment. Pin the dependency 'mcp&lt;2': the upstream package leaves mcp unpinned and its server fails to import under mcp 2.x (FastMCP was renamed). |
| 副作用 | Read-only HTTPS queries to NCBI E-utilities (Gene, ClinVar, Taxonomy, Nucleotide, PubMed) and the UniProt REST API; every request sends tool and a fixed contact email parameter to NCBI, nothing is written locally. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`get_gene_details`](#get_gene_details) | `genome_annotation` | Retrieves detailed NCBI Gene information (name, summary, location) for an NCBI Gene ID. |
| [`get_ncbi_sequence`](#get_ncbi_sequence) | `sequence` | Fetches a nucleotide or protein sequence from NCBI by accession in FASTA or GenBank format. |
| [`search_clinvar_variants`](#search_clinvar_variants) | `variant` | Searches NCBI ClinVar for variants in a gene filtered by clinical significance (pathogenic, benign, uncertain, all). |
| [`search_genes`](#search_genes) | `genome_annotation` | Searches the NCBI Gene database by symbol or keyword for a given organism. |

## `get_gene_details` {#get_gene_details}

Retrieves detailed NCBI Gene information (name, summary, location) for an NCBI Gene ID.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`gene_id`
- **关键词：**ncbi gene, entrez, gene id, gene summary, genomic location

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 |  |

??? note "上游说明（英文原文）"

    Ruft detaillierte Informationen zu einem Gen via NCBI Gene ID ab.

    Args:
    gene_id: NCBI Gene ID (z.B. "672" für BRCA1, "7157" für TP53)

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-ncbi-uniprot',
    'get_gene_details',
    {'gene_id': '...'},
)
```

## `get_ncbi_sequence` {#get_ncbi_sequence}

Fetches a nucleotide or protein sequence from NCBI by accession in FASTA or GenBank format.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`accession`
- **关键词：**ncbi, accession, fasta, genbank, nucleotide, mrna

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 |  |
| `database` | string | 否 | 默认 `nucleotide` |
| `sequence_type` | string | 否 | 默认 `fasta` |

??? note "上游说明（英文原文）"

    Ruft eine DNA- oder Protein-Sequenz aus NCBI ab.

    Args:
    accession: NCBI Accession-Nummer (z.B. "NM_007294" für BRCA1 mRNA, "NP_009225" für BRCA1 Protein)
    database: Datenbank — "nucleotide" (DNA/RNA) oder "protein". Standard: nucleotide
    sequence_type: Format — "fasta" (Sequenz) oder "gb" (GenBank-Eintrag). Standard: fasta

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-ncbi-uniprot',
    'get_ncbi_sequence',
    {'accession': '...'},
)
```

## `search_clinvar_variants` {#search_clinvar_variants}

Searches NCBI ClinVar for variants in a gene filtered by clinical significance (pathogenic, benign, uncertain, all).

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene_symbol`
- **关键词：**clinvar, pathogenic, clinical significance, mutation, gene variants

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 |  |
| `significance` | string | 否 | 默认 `pathogenic` |
| `max_results` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Sucht klinische Varianten (Mutationen) für ein Gen via NCBI ClinVar.

    Args:
    gene_symbol: Gen-Symbol (z.B. "BRCA1", "TP53", "CFTR")
    significance: Klinische Signifikanz — "pathogenic", "benign", "uncertain" oder "all". Standard: pathogenic
    max_results: Maximale Ergebnisanzahl (1-20). Standard: 10

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-ncbi-uniprot',
    'search_clinvar_variants',
    {'gene_symbol': '...'},
)
```

## `search_genes` {#search_genes}

Searches the NCBI Gene database by symbol or keyword for a given organism.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`query`
- **关键词：**ncbi gene, gene search, gene symbol, entrez, organism

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 |  |
| `organism` | string | 否 | 默认 `human` |
| `max_results` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Sucht nach Genen in der NCBI Gene Datenbank.

    Args:
    query: Suchbegriff (z.B. "BRCA1", "TP53", "insulin receptor")
    organism: Organismus (z.B. "human", "mouse", "rat"). Standard: human
    max_results: Maximale Ergebnisanzahl (1-20). Standard: 10

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-ncbi-uniprot',
    'search_genes',
    {'query': '...'},
)
```
