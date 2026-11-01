# NCBI Database MCP (hpend2373) (`ncbi-database`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **4** 个 DNA/基因相关工具（服务器共列出 4 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [hpend2373/NCBI-Database-MCP @ `03a6b80`](https://github.com/hpend2373/NCBI-Database-MCP/tree/03a6b80490753ba47166903daa7d7bce35848af3) |
| 许可证 | MIT declared in pyproject.toml; standalone license text not found |
| 启动命令 | `python src/gene_to_genomic_server.py` |
| 传输方式 | stdio |
| 安装 | Prepare the pinned NCBI-Database-MCP checkout and its requirements, set workdir to the checkout and executable_path to its Python. Uses the README-supported standard MCP script rather than executing its shell installer/launcher. |
| 副作用 | Queries NCBI gene/sequence/GEO resources; selected results can contain supplied gene or region queries. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`gene_to_genomic_sequence`](#gene_to_genomic_sequence) | `sequence` | Retrieves the genomic, CDS, mRNA or protein sequence of a gene by symbol and organism via NCBI E-utilities. |
| [`get_genomic_sequence`](#get_genomic_sequence) | `sequence` | Retrieves genomic DNA sequence for a chromosome accession (e.g. NC_000017.11) between start and end positions. |
| [`search_gene_info`](#search_gene_info) | `genome_annotation` | Searches NCBI Gene for a gene by name and organism and returns its information and genomic location. |
| [`search_geo_datasets`](#search_geo_datasets) | `sequencing_data` | Searches NCBI GEO datasets by disease and organism, filtered by study type such as RNA-seq, array or ChIP-seq. |

## `gene_to_genomic_sequence` {#gene_to_genomic_sequence}

Retrieves the genomic, CDS, mRNA or protein sequence of a gene by symbol and organism via NCBI E-utilities.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`gene_name`
- **关键词：**gene symbol, cds, mrna, fasta, genbank

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_name` | string | 是 | Gene symbol or name (e.g., TP53, BRCA1) |
| `organism` | string | 否 | Organism name (default: human) 默认 `human` |
| `sequence_type` | string | 否 | Type of sequence to retrieve (default: genomic) 可选值：`genomic`, `cds`, `mrna`, `protein` 默认 `genomic` |
| `output_format` | string | 否 | Output format (default: fasta) 可选值：`fasta`, `genbank`, `json` 默认 `fasta` |

??? note "上游说明（英文原文）"

    Convert gene name to genomic DNA sequence using NCBI E-utilities

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-database',
    'gene_to_genomic_sequence',
    {'gene_name': '...'},
)
```

## `get_genomic_sequence` {#get_genomic_sequence}

Retrieves genomic DNA sequence for a chromosome accession (e.g. NC_000017.11) between start and end positions.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`chromosome`, `start`, `end`
- **关键词：**chromosome, coordinates, refseq, fasta, genbank

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome accession (e.g., NC_000017.11) |
| `start` | integer | 是 | Start position |
| `end` | integer | 是 | End position |
| `output_format` | string | 否 | Output format (default: fasta) 可选值：`fasta`, `genbank` 默认 `fasta` |

??? note "上游说明（英文原文）"

    Get genomic sequence from chromosome coordinates

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-database',
    'get_genomic_sequence',
    {'chromosome': '...', 'start': 0, 'end': 0},
)
```

## `search_gene_info` {#search_gene_info}

Searches NCBI Gene for a gene by name and organism and returns its information and genomic location.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`gene_name`
- **关键词：**ncbi gene, genomic location, gene id, chromosome, organism

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_name` | string | 是 | Gene symbol or name |
| `organism` | string | 否 | Organism name (default: human) 默认 `human` |

??? note "上游说明（英文原文）"

    Search for gene information including genomic location

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-database',
    'search_gene_info',
    {'gene_name': '...'},
)
```

## `search_geo_datasets` {#search_geo_datasets}

Searches NCBI GEO datasets by disease and organism, filtered by study type such as RNA-seq, array or ChIP-seq.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`disease`
- **关键词：**gse, expression profiling, rna-seq, microarray, disease

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `disease` | string | 是 | Disease or condition name (e.g., 'cancer', 'diabetes', 'Alzheimer') |
| `organism` | string | 否 | Organism (default: 'Homo sapiens') 可选值：`Homo sapiens`, `Mus musculus`, `Rattus norvegicus` 默认 `Homo sapiens` |
| `study_type` | string | 否 | Type of expression study (default: Expression profiling by high throughput sequencing) 可选值：`Expression profiling by array`, `Expression profiling by high throughput sequencing`, `Genome binding/occupancy profiling by high throughput sequencing` 默认 `Expression profiling by high throughput sequencing` |
| `max_results` | integer | 否 | Maximum number of results to return (default: 10) 默认 `10` |

??? note "上游说明（英文原文）"

    Search GEO datasets by disease/condition and organism

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'ncbi-database',
    'search_geo_datasets',
    {'disease': '...'},
)
```
