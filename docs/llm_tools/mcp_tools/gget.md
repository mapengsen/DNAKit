# gget MCP (`gget`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **15** 个 DNA/基因相关工具（服务器共列出 20 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [longevity-genie/gget-mcp @ `f3c6312`](https://github.com/longevity-genie/gget-mcp/tree/f3c63125e29954c20b3a81d1d5df70405ddf5896) |
| 许可证 | MIT |
| 启动命令 | `gget-mcp stdio` |
| 传输方式 | stdio |
| 安装 | Install gget-mcp and its separately specified gget dependencies in Python &gt;=3.10; the verified command is gget-mcp stdio. |
| 副作用 | Queries gget-supported services; selected tools may run external searches, load models, download assets or write outputs. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`gget_archs4`](#gget_archs4) | `expression` | Finds the most correlated genes or the tissue expression atlas for a gene using ARCHS4 RNA-seq data. |
| [`gget_bgee`](#gget_bgee) | `expression` | Retrieves orthologs or anatomical expression data for an Ensembl gene ID from the Bgee database. |
| [`gget_blast`](#gget_blast) | `sequence` | Runs NCBI BLAST (blastn, blastp, blastx, tblastn, tblastx) for a sequence against nt, nr, RefSeq or SwissProt. |
| [`gget_blat`](#gget_blat) | `sequence` | Aligns a nucleotide or amino acid sequence to a UCSC genome assembly (hg38, mm39, taeGut2) with BLAT to get genomic location. |
| [`gget_cellxgene`](#gget_cellxgene) | `expression` | Queries CZ CELLxGENE Census single-cell RNA-seq data filtered by genes, tissues and cell types. |
| [`gget_cosmic`](#gget_cosmic) | `variant` | Searches the COSMIC database for cancer somatic mutations by gene symbol or name. |
| [`gget_enrichr`](#gget_enrichr) | `genome_annotation` | Performs Enrichr enrichment analysis on a gene list against KEGG, ChEA, GO, GWAS, PanglaoDB or KEA libraries. |
| [`gget_info`](#gget_info) | `genome_annotation` | Retrieves gene and transcript metadata (from Ensembl, UniProt, NCBI) for one or more Ensembl IDs. |
| [`gget_muscle`](#gget_muscle) | `sequence` | Aligns multiple nucleotide or amino acid sequences with MUSCLE and saves the result as FASTA/aligned FASTA. |
| [`gget_mutate`](#gget_mutate) | `sequence` | Mutates nucleotide sequences according to mutations in standard annotation (e.g. c.2C&gt;T). |
| [`gget_opentargets`](#gget_opentargets) | `association` | Queries Open Targets for a gene's associated diseases, drugs, tractability, pharmacogenetics, expression, DepMap or interactions. |
| [`gget_ref`](#gget_ref) | `genome_annotation` | Fetches FTP download links for Ensembl reference genome DNA, cDNA, CDS, GTF or peptide files for a species. |
| [`gget_search`](#gget_search) | `genome_annotation` | Searches Ensembl for genes or transcripts whose names or descriptions match free-text terms in a species. |
| [`gget_search_genes`](#gget_search_genes) | `genome_annotation` | Searches Ensembl by gene symbol(s) with descriptive terms to return Ensembl gene or transcript IDs. |
| [`gget_seq`](#gget_seq) | `sequence` | Fetches nucleotide sequences (or translated amino acid sequences) for Ensembl IDs and saves them as FASTA. |

## `gget_archs4` {#gget_archs4}

Finds the most correlated genes or the tissue expression atlas for a gene using ARCHS4 RNA-seq data.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene`
- **关键词：**rna-seq, coexpression, tissue expression, correlated genes

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 | Gene symbol (e.g., 'STAT4') or Ensembl ID if ensembl=True (e.g., 'ENSG00000138378') |
| `which` | string | 否 | Analysis type - 'correlation' (most correlated genes) or 'tissue' (tissue expression atlas) 默认 `correlation` |
| `species` | string | 否 | Target species - 'human' or 'mouse' (only for tissue expression atlas) 默认 `human` |

??? note "上游说明（英文原文）"

    Find correlated genes or tissue expression atlas using ARCHS4 RNA-seq database.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_archs4',
    {'gene': '...'},
)
```

## `gget_bgee` {#gget_bgee}

Retrieves orthologs or anatomical expression data for an Ensembl gene ID from the Bgee database.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`gene_id`
- **关键词：**bgee, orthologs, tissue expression, ensembl id, anatomy

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 | Ensembl gene ID (e.g., 'ENSG00000012048' for BRCA1) |
| `type` | string | 否 | Type of data to retrieve - 'orthologs' or 'expression' 默认 `orthologs` |

??? note "上游说明（英文原文）"

    Get orthologs or expression data for a gene from Bgee database.

    PREREQUISITE: Use search_genes to get Ensembl ID first.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_bgee',
    {'gene_id': '...'},
)
```

## `gget_blast` {#gget_blast}

Runs NCBI BLAST (blastn, blastp, blastx, tblastn, tblastx) for a sequence against nt, nr, RefSeq or SwissProt.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**blastn, similarity search, nt, refseq, alignment

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | Nucleotide or amino acid sequence (string) or path to FASTA file |
| `program` | string | 否 | BLAST program - 'blastn', 'blastp', 'blastx', 'tblastn', 'tblastx', or 'default' (auto-detect) 默认 `default` |
| `database` | string | 否 | BLAST database - 'nt', 'nr', 'refseq_rna', 'refseq_protein', 'swissprot', or 'default' (auto-detect) 默认 `default` |

??? note "上游说明（英文原文）"

    BLAST a nucleotide or amino acid sequence against any BLAST database.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_blast',
    {'sequence': '...'},
)
```

## `gget_blat` {#gget_blat}

Aligns a nucleotide or amino acid sequence to a UCSC genome assembly (hg38, mm39, taeGut2) with BLAT to get genomic location.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**ucsc, genome mapping, hg38, mm39, genomic coordinates

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | Nucleotide or amino acid sequence (string) or path to FASTA file containing one sequence |
| `assembly` | string | 否 | Genome assembly - 'human' (hg38), 'mouse' (mm39), 'zebrafinch' (taeGut2) 默认 `human` |

??? note "上游说明（英文原文）"

    BLAT a nucleotide or amino acid sequence against any BLAT UCSC assembly.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_blat',
    {'sequence': '...'},
)
```

## `gget_cellxgene` {#gget_cellxgene}

Queries CZ CELLxGENE Census single-cell RNA-seq data filtered by genes, tissues and cell types.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**无
- **关键词：**single-cell, scrna-seq, census, cell type, anndata

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | array \| string \| null | 否 | Gene name(s) or Ensembl ID(s) (e.g., \['ACE2', 'SLC5A1'\]) 默认 `None` |
| `tissue` | array \| string \| null | 否 | Tissue(s) to query (e.g., \['lung', 'blood'\]) 默认 `None` |
| `cell_type` | array \| string \| null | 否 | Cell type(s) to query (e.g., \['mucus secreting cell'\]) 默认 `None` |
| `species` | string | 否 | Target species - 'homo_sapiens' or 'mus_musculus' 默认 `homo_sapiens` |

??? note "上游说明（英文原文）"

    Query single-cell RNA-seq data from CZ CELLxGENE Discover using Census.

    NOTE: Querying large datasets requires &gt;16 GB RAM and &gt;5 Mbps internet connection.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_cellxgene',
    {},
)
```

## `gget_cosmic` {#gget_cosmic}

Searches the COSMIC database for cancer somatic mutations by gene symbol or name.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`searchterm`
- **关键词：**cosmic, somatic mutation, cancer, mutations

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `searchterm` | string | 是 | Gene symbol or name to search for (e.g., 'PIK3CA', 'BRCA1') |
| `cosmic_tsv_path` | string \| null | 否 | Path to COSMIC TSV file (optional, uses default if None) 默认 `None` |

??? note "上游说明（英文原文）"

    Search COSMIC database for cancer mutations or download COSMIC databases.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_cosmic',
    {'searchterm': '...'},
)
```

## `gget_enrichr` {#gget_enrichr}

Performs Enrichr enrichment analysis on a gene list against KEGG, ChEA, GO, GWAS, PanglaoDB or KEA libraries.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`genes`
- **关键词：**enrichment, gene set, kegg, go, pathway

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genes` | array of string | 是 | List of gene symbols (e.g., \['PHF14', 'RBM3'\]) or Ensembl IDs if ensembl=True |
| `database` | string | 否 | Reference database shortcuts: 'pathway' (KEGG), 'transcription' (ChEA), 'ontology' (GO), 'diseases_drugs' (GWAS), 'celltypes' (PanglaoDB), 'kinase_interactions' (KEA) 默认 `pathway` |
| `species` | string | 否 | Species database - 'human', 'mouse', 'fly', 'yeast', 'worm', 'fish' 默认 `human` |

??? note "上游说明（英文原文）"

    Perform functional enrichment analysis on gene list using Enrichr.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_enrichr',
    {'genes': []},
)
```

## `gget_info` {#gget_info}

Retrieves gene and transcript metadata (from Ensembl, UniProt, NCBI) for one or more Ensembl IDs.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`ensembl_ids`
- **关键词：**ensembl id, transcript, gene metadata, uniprot, ncbi

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_ids` | string \| array | 是 | One or more Ensembl gene IDs (e.g., 'ENSG00000141510' or \['ENSG00000141510'\]) Also supports WormBase and FlyBase IDs |

??? note "上游说明（英文原文）"

    Get detailed gene and transcript metadata using Ensembl IDs.

    PREREQUISITE: Use search_genes first to get Ensembl IDs from gene names/symbols.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_info',
    {'ensembl_ids': '...'},
)
```

## `gget_muscle` {#gget_muscle}

Aligns multiple nucleotide or amino acid sequences with MUSCLE and saves the result as FASTA/aligned FASTA.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequences`
- **关键词：**multiple sequence alignment, msa, fasta, afa

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequences` | array \| string | 是 | List of sequences or ABSOLUTE path to FASTA file containing sequences to be aligned |
| `output_path` | string \| null | 否 | ABSOLUTE path to output file (e.g., '/home/user/alignment.fasta'). AVOID relative paths as they cause file location issues. Auto-generated if not provided. 默认 `None` |
| `format` | string | 否 | Output format ('fasta' for FASTA format, 'afa' for aligned FASTA format) 可选值：`fasta`, `afa` 默认 `fasta` |

??? note "上游说明（英文原文）"

    Align sequences using MUSCLE and save to local file.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_muscle',
    {'sequences': '...'},
)
```

## `gget_mutate` {#gget_mutate}

Mutates nucleotide sequences according to mutations in standard annotation (e.g. c.2C&gt;T).

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequences`, `mutations`
- **关键词：**mutation, hgvs, nucleotide sequence, substitution, in silico mutagenesis

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequences` | string \| array | 是 | Path to FASTA file or sequences as string/list (e.g., 'AGCTAGCT' or \['ACTG', 'ATCG'\]) |
| `mutations` | string \| array | 是 | Path to CSV/TSV file or mutations as string/list (e.g., 'c.2C&gt;T' or \['c.2C&gt;T', 'c.1A&gt;C'\]) Standard mutation annotation format required |

??? note "上游说明（英文原文）"

    Mutate nucleotide sequences according to provided mutations in standard annotation.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_mutate',
    {'sequences': '...', 'mutations': '...'},
)
```

## `gget_opentargets` {#gget_opentargets}

Queries Open Targets for a gene's associated diseases, drugs, tractability, pharmacogenetics, expression, DepMap or interactions.

- **类别：**`association`（GWAS, QTL, LD and gene/variant-disease or trait associations）
- **必填参数：**`ensembl_id`
- **关键词：**open targets, gene-disease, target, depmap, tractability

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 | Ensembl gene ID (e.g., 'ENSG00000169194') |
| `resource` | string | 否 | Type of information - 'diseases', 'drugs', 'tractability', 'pharmacogenetics', 'expression', 'depmap', 'interactions' 默认 `diseases` |

??? note "上游说明（英文原文）"

    Query OpenTargets for diseases, drugs, and other data associated with a gene.

    PREREQUISITE: Use search_genes to get Ensembl ID first.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_opentargets',
    {'ensembl_id': '...'},
)
```

## `gget_ref` {#gget_ref}

Fetches FTP download links for Ensembl reference genome DNA, cDNA, CDS, GTF or peptide files for a species.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**无
- **关键词：**reference genome, gtf, ftp, cdna, ensembl release

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `species` | string | 否 | Species in format "genus_species" (e.g., "homo_sapiens"). Shortcuts supported: "human", "mouse", "human_grch37" 默认 `homo_sapiens` |
| `which` | string \| array | 否 | Which results to return. Options: 'gtf', 'cdna', 'dna', 'cds', 'cdrna', 'pep', 'all' 默认 `all` |

??? note "上游说明（英文原文）"

    Fetch FTPs for reference genomes and annotations by species from Ensembl.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_ref',
    {},
)
```

## `gget_search` {#gget_search}

Searches Ensembl for genes or transcripts whose names or descriptions match free-text terms in a species.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`search_terms`
- **关键词：**ensembl, gene search, free text, synonyms, species

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `search_terms` | string \| array | 是 | Search terms, names, or synonyms (e.g., 'cancer' or \['apoptosis', 'death'\]) |
| `species` | string | 否 | Target species (e.g., 'homo_sapiens', 'mus_musculus') 默认 `homo_sapiens` |

??? note "上游说明（英文原文）"

    General search for any biological terms using gene symbols, names, or synonyms.

    This is a general search that looks broadly across gene names and descriptions.
    For specific gene symbol searches, use search_genes_simple instead.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_search',
    {'search_terms': '...'},
)
```

## `gget_search_genes` {#gget_search_genes}

Searches Ensembl by gene symbol(s) with descriptive terms to return Ensembl gene or transcript IDs.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`search_terms`
- **关键词：**ensembl id, gene symbol, id resolution, batch

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `search_terms` | string \| array | 是 | SINGLE gene symbol OR LIST of gene symbols with descriptive terms Single: 'APP amyloid precursor' Batch: \['BACE1 beta secretase', 'MAPT tau', 'APOE apolipoprotein'\] |
| `species` | string | 否 | Target species (e.g., 'homo_sapiens', 'mus_musculus') 默认 `homo_sapiens` |
| `id_type` | string | 否 | "gene" (default) or "transcript" - whether to return genes or transcripts 默认 `gene` |

??? note "上游说明（英文原文）"

    Search for specific genes using gene symbols with enhanced search strategy.

    🚀 **BATCH PROCESSING SUPPORTED**: This function can process multiple genes in a single call!
    Use this tool FIRST when you have gene names/symbols and need to find their Ensembl IDs.
    Returns Ensembl IDs which are required for get_gene_info and get_sequences tools.

    IMPORTANT: Due to limitations in Ensembl search, short gene names often fail to find results.
    For best results, provide descriptive terms along with gene symbols:

    RECOMMENDED FORMAT: "GENE_SYMBOL descriptive_terms"
    Examples:
    - Instead of: "APP"
    - Use: "APP amyloid precursor" or "APP amyloid beta precursor protein"
    - Instead of: \["BACE1", "MAPT"\]
    - Use: \["BACE1 beta secretase", "MAPT microtubule tau"\]

    This function uses AND search for multi-word terms and OR search for single words.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_search_genes',
    {'search_terms': '...'},
)
```

## `gget_seq` {#gget_seq}

Fetches nucleotide sequences (or translated amino acid sequences) for Ensembl IDs and saves them as FASTA.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`ensembl_ids`
- **关键词：**fasta, ensembl id, nucleotide sequence, transcript sequence

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_ids` | string \| array | 是 | One or more Ensembl gene IDs (e.g., 'ENSG00000141510' or \['ENSG00000141510'\]) |
| `translate` | boolean | 否 | If True, returns amino acid sequences; if False, returns nucleotide sequences 默认 `False` |
| `output_path` | string \| null | 否 | ABSOLUTE path to output file (e.g., '/home/user/sequences.fasta'). AVOID relative paths as they cause file location issues. Auto-generated if not provided. 默认 `None` |
| `format` | string | 否 | Output format (currently supports 'fasta') 默认 `fasta` |

??? note "上游说明（英文原文）"

    Fetch sequences and save to local file in stdio mode.

    PREREQUISITE: Use search_genes first to get Ensembl IDs from gene names/symbols.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'gget',
    'gget_seq',
    {'ensembl_ids': '...'},
)
```
