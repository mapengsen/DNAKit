# AWS Kiro Life Sciences: genomics (`kiro-genomics`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **17** 个 DNA/基因相关工具（服务器共列出 18 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [aws-samples/sample-kiro-power-life-sciences @ `95b67e7`](https://github.com/aws-samples/sample-kiro-power-life-sciences/tree/95b67e763df43e94e5bcfd8a230a58553a888cfc) |
| 许可证 | MIT-0 in root LICENSE.txt; MIT declared by component pyproject.toml |
| 启动命令 | `life-sciences-genomics` |
| 传输方式 | stdio |
| 安装 | Prepare the pinned aws-samples/sample-kiro-power-life-sciences source. Separately install life-sciences-common and the life-sciences-genomics component in its Python &gt;=3.10 environment, then set executable_path to the installed life-sciences-genomics console script. No uvx or automatic install is used. |
| 副作用 | Selected tools query external biological resources; molecular-analysis and workflow tools can submit sequences or access local result files. Named credentials are forwarded only for providers that need them. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`clinvar_get_variation`](#clinvar_get_variation) | `variant` | Fetches the full ClinVar variation record for a given ClinVar variation ID. |
| [`clinvar_search`](#clinvar_search) | `variant` | Searches ClinVar by rsID, HGVS notation or gene name and returns matching variation records. |
| [`cosmic_search`](#cosmic_search) | `variant` | Searches the COSMIC database for somatic cancer mutations in a given gene symbol (requires COSMIC_API_KEY). |
| [`dbsnp_lookup`](#dbsnp_lookup) | `variant` | Retrieves the dbSNP record for a variant given its rsID. |
| [`ddbj_fetch_sequence`](#ddbj_fetch_sequence) | `sequence` | Fetches a sequence record from DDBJ by accession number in FASTA format. |
| [`ddbj_search`](#ddbj_search) | `sequence` | Searches DDBJ (DNA Data Bank of Japan) sequence entries by keyword. |
| [`encode_search`](#encode_search) | `sequencing_data` | Searches ENCODE experiments filtered by biosample, assay and/or target label. |
| [`ensembl_gene_lookup`](#ensembl_gene_lookup) | `genome_annotation` | Looks up an Ensembl gene record (ID, location, biotype) by gene symbol and species via the Ensembl REST API. |
| [`ensembl_sequence`](#ensembl_sequence) | `sequence` | Retrieves the nucleotide sequence for a genomic region of a species via the Ensembl REST API. |
| [`ensembl_variants`](#ensembl_variants) | `variant` | Lists known variants overlapping a genomic region (e.g. 7:140424943-140624564) for a species via Ensembl REST. |
| [`geo_get_dataset`](#geo_get_dataset) | `sequencing_data` | Retrieves a GEO dataset/series record by accession (e.g. GSE12345). |
| [`geo_search`](#geo_search) | `sequencing_data` | Searches NCBI GEO datasets by keyword and returns matching dataset summaries. |
| [`gnomad_variant`](#gnomad_variant) | `variant` | Retrieves gnomAD allele frequency data for a variant given in chrom-pos-ref-alt form. |
| [`ncbi_fetch_sequence`](#ncbi_fetch_sequence) | `sequence` | Fetches a nucleotide sequence record by accession via NCBI Entrez efetch in FASTA or GenBank format. |
| [`ncbi_search`](#ncbi_search) | `genome_annotation` | Runs an NCBI Entrez esearch over a chosen database (e.g. gene, nucleotide, protein) and returns matching IDs. |
| [`sra_search`](#sra_search) | `sequencing_data` | Searches NCBI SRA runs by keyword or BioProject accession. |
| [`thousand_genomes_frequency`](#thousand_genomes_frequency) | `variant` | Retrieves 1000 Genomes population allele frequencies for a variant (e.g. rsID) via Ensembl. |

## `clinvar_get_variation` {#clinvar_get_variation}

Fetches the full ClinVar variation record for a given ClinVar variation ID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variation_id`
- **关键词：**clinical significance, pathogenicity, variation id, germline variant, hgvs

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variation_id` | string | 是 | ClinVar variation ID. |

??? note "上游说明（英文原文）"

    Get a full ClinVar variation record by variation ID.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'clinvar_get_variation',
    {'variation_id': '...'},
)
```

## `clinvar_search` {#clinvar_search}

Searches ClinVar by rsID, HGVS notation or gene name and returns matching variation records.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`query`
- **关键词：**rsid, hgvs, pathogenic, clinical significance, gene variants

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search term (rsID, HGVS, or gene name). |
| `max_results` | integer | 否 | Maximum results. 默认 `10` |

??? note "上游说明（英文原文）"

    Search ClinVar by rsID, HGVS notation, or gene name.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'clinvar_search',
    {'query': '...'},
)
```

## `cosmic_search` {#cosmic_search}

Searches the COSMIC database for somatic cancer mutations in a given gene symbol (requires COSMIC_API_KEY).

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene`
- **关键词：**somatic mutation, cancer, tumor, mutation catalogue, oncogene

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (e.g. 'BRAF'). |

??? note "上游说明（英文原文）"

    Search COSMIC somatic mutations by gene name. Requires COSMIC_API_KEY.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'cosmic_search',
    {'gene': '...'},
)
```

## `dbsnp_lookup` {#dbsnp_lookup}

Retrieves the dbSNP record for a variant given its rsID.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`rsid`
- **关键词：**rsid, snp, single nucleotide polymorphism, refsnp, allele

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `rsid` | string | 是 | dbSNP rsID (e.g. 'rs328'). |

??? note "上游说明（英文原文）"

    Lookup a variant by rsID from dbSNP.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'dbsnp_lookup',
    {'rsid': '...'},
)
```

## `ddbj_fetch_sequence` {#ddbj_fetch_sequence}

Fetches a sequence record from DDBJ by accession number in FASTA format.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`accession`
- **关键词：**fasta, nucleotide, accession, insdc, dna sequence

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | DDBJ accession number. |

??? note "上游说明（英文原文）"

    Fetch a sequence from DDBJ in FASTA format.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ddbj_fetch_sequence',
    {'accession': '...'},
)
```

## `ddbj_search` {#ddbj_search}

Searches DDBJ (DNA Data Bank of Japan) sequence entries by keyword.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`query`
- **关键词：**nucleotide, insdc, accession, sequence records, dna databank

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 | Search keyword. |
| `max_results` | integer | 否 | Maximum results. 默认 `10` |

??? note "上游说明（英文原文）"

    Search DDBJ by keyword.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ddbj_search',
    {'query': '...'},
)
```

## `encode_search` {#encode_search}

Searches ENCODE experiments filtered by biosample, assay and/or target label.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**无
- **关键词：**chip-seq, atac-seq, biosample, assay, functional genomics

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `biosample` | string | 否 | Biosample term name. 默认 `` |
| `assay` | string | 否 | Assay title. 默认 `` |
| `target` | string | 否 | Target label. 默认 `` |

??? note "上游说明（英文原文）"

    Search ENCODE experiments by biosample, assay, and/or target.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'encode_search',
    {},
)
```

## `ensembl_gene_lookup` {#ensembl_gene_lookup}

Looks up an Ensembl gene record (ID, location, biotype) by gene symbol and species via the Ensembl REST API.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`symbol`
- **关键词：**gene symbol, ensembl id, genomic location, biotype, species

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `symbol` | string | 是 | Gene symbol (e.g. 'BRCA2'). |
| `species` | string | 否 | Species name. 默认 `homo_sapiens` |

??? note "上游说明（英文原文）"

    Lookup a gene by symbol and species via Ensembl REST.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ensembl_gene_lookup',
    {'symbol': '...'},
)
```

## `ensembl_sequence` {#ensembl_sequence}

Retrieves the nucleotide sequence for a genomic region of a species via the Ensembl REST API.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`region`
- **关键词：**genomic region, dna sequence, reference genome, coordinates, fasta

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `region` | string | 是 | Genomic region. |
| `species` | string | 否 | Species name. 默认 `homo_sapiens` |

??? note "上游说明（英文原文）"

    Get nucleotide sequence for a genomic region via Ensembl REST.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ensembl_sequence',
    {'region': '...'},
)
```

## `ensembl_variants` {#ensembl_variants}

Lists known variants overlapping a genomic region (e.g. 7:140424943-140624564) for a species via Ensembl REST.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`region`
- **关键词：**genomic region, snp, overlap, variation, coordinates

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `region` | string | 是 | Genomic region (e.g. '7:140424943-140624564'). |
| `species` | string | 否 | Species name. 默认 `homo_sapiens` |

??? note "上游说明（英文原文）"

    Get variants in a genomic region via Ensembl REST.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ensembl_variants',
    {'region': '...'},
)
```

## `geo_get_dataset` {#geo_get_dataset}

Retrieves a GEO dataset/series record by accession (e.g. GSE12345).

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`accession`
- **关键词：**gse, expression dataset, series, microarray, rna-seq

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | GEO accession number. |

??? note "上游说明（英文原文）"

    Get a GEO dataset by accession (e.g. GSE12345).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'geo_get_dataset',
    {'accession': '...'},
)
```

## `geo_search` {#geo_search}

Searches NCBI GEO datasets by keyword and returns matching dataset summaries.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`term`
- **关键词：**gse, expression data, rna-seq, microarray, datasets

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `term` | string | 是 | Search keyword. |
| `max_results` | integer | 否 | Maximum results. 默认 `10` |

??? note "上游说明（英文原文）"

    Search GEO datasets by keyword.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'geo_search',
    {'term': '...'},
)
```

## `gnomad_variant` {#gnomad_variant}

Retrieves gnomAD allele frequency data for a variant given in chrom-pos-ref-alt form.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant_id`
- **关键词：**allele frequency, population frequency, gnomad, chrom-pos-ref-alt, maf

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant_id` | string | 是 | Variant identifier (e.g. '1-55516888-G-A'). |

??? note "上游说明（英文原文）"

    Get allele frequency data for a variant from gnomAD.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'gnomad_variant',
    {'variant_id': '...'},
)
```

## `ncbi_fetch_sequence` {#ncbi_fetch_sequence}

Fetches a nucleotide sequence record by accession via NCBI Entrez efetch in FASTA or GenBank format.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`accession`
- **关键词：**efetch, genbank, fasta, nucleotide, accession

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 | NCBI accession number. |
| `format` | string | 否 | Output format ('fasta' or 'gb'). 默认 `fasta` |

??? note "上游说明（英文原文）"

    Fetch a nucleotide sequence by accession via Entrez efetch.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ncbi_fetch_sequence',
    {'accession': '...'},
)
```

## `ncbi_search` {#ncbi_search}

Runs an NCBI Entrez esearch over a chosen database (e.g. gene, nucleotide, protein) and returns matching IDs.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`database`, `term`
- **关键词：**entrez, esearch, gene database, nucleotide, uid

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `database` | string | 是 | NCBI database name (e.g. 'gene', 'nucleotide', 'protein'). |
| `term` | string | 是 | Search query string. |
| `max_results` | integer | 否 | Maximum results to return. 默认 `20` |

??? note "上游说明（英文原文）"

    Search any NCBI database via Entrez esearch.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'ncbi_search',
    {'database': '...', 'term': '...'},
)
```

## `sra_search` {#sra_search}

Searches NCBI SRA runs by keyword or BioProject accession.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`term`
- **关键词：**sequencing runs, bioproject, reads, srr, raw data

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `term` | string | 是 | Search keyword or BioProject accession. |
| `max_results` | integer | 否 | Maximum results. 默认 `10` |

??? note "上游说明（英文原文）"

    Search SRA runs by keyword or BioProject.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'sra_search',
    {'term': '...'},
)
```

## `thousand_genomes_frequency` {#thousand_genomes_frequency}

Retrieves 1000 Genomes population allele frequencies for a variant (e.g. rsID) via Ensembl.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`variant`
- **关键词：**1000 genomes, allele frequency, population, rsid, superpopulation

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant` | string | 是 | Variant identifier (e.g. 'rs56116432'). |
| `species` | string | 否 | Species name. 默认 `homo_sapiens` |

??? note "上游说明（英文原文）"

    Get population allele frequencies from 1000 Genomes.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'kiro-genomics',
    'thousand_genomes_frequency',
    {'variant': '...'},
)
```
