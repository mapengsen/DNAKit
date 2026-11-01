# bioinformatics-mcp (`bioinformatics-mcp`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **33** 个 DNA/基因相关工具（服务器共列出 70 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [raonyguimaraes/bioinformatics-mcp @ `b607d56`](https://github.com/raonyguimaraes/bioinformatics-mcp/tree/b607d569a95fad247125894efa4b3e14e024c5e2) |
| 许可证 | MIT |
| 启动命令 | `bioinformatics-mcp` |
| 传输方式 | stdio |
| 安装 | Install separately into a Python &gt;=3.11 environment from the pinned source (pip install "bioinformatics-mcp\[full\] @ git+https://github.com/raonyguimaraes/bioinformatics-mcp.git@b607d569a95fad247125894efa4b3e14e024c5e2"; the README's PyPI name returned 404 on PyPI); no environment variables are used. |
| 副作用 | Database tools send queries to NCBI E-utilities, UniProt, Ensembl and STRING over HTTPS; save_anndata and export_obs_table write files to caller-supplied paths; other tools read local sequence, VCF, BAM, CSV and AnnData files and hold datasets in memory. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`bam_stats`](#bam_stats) | `sequencing_data` | Computes summary statistics (mapped/unmapped reads, references) from a BAM or SAM alignment file. |
| [`cpg_density`](#cpg_density) | `sequence` | Computes CpG observed/expected ratio in sliding windows along a DNA sequence and classifies windows as island, shore or open sea. |
| [`export_gene_expression`](#export_gene_expression) | `expression` | Extracts expression values for a list of genes across all cells of an AnnData dataset, optionally from a layer. |
| [`fetch_ensembl_sequence`](#fetch_ensembl_sequence) | `sequence` | Fetches a genomic, cDNA or CDS nucleotide sequence from Ensembl for a gene or transcript ID. |
| [`fetch_gene_ontology`](#fetch_gene_ontology) | `genome_annotation` | Fetches GO biological process, molecular function and cellular component annotations for a gene symbol via UniProt. |
| [`fetch_ncbi_gene_details`](#fetch_ncbi_gene_details) | `genome_annotation` | Fetches NCBI Gene details (name, summary, genomic location, references) for an NCBI Gene ID. |
| [`fetch_ncbi_sequence`](#fetch_ncbi_sequence) | `sequence` | Downloads a nucleotide sequence from NCBI by accession in FASTA or GenBank format, with optional truncation. |
| [`find_differentially_methylated_regions`](#find_differentially_methylated_regions) | `regulation_epigenome` | Identifies differentially methylated regions between two sample groups using a sliding genomic window and a minimum methylation difference. |
| [`find_orfs`](#find_orfs) | `sequence` | Finds all open reading frames above a minimum length in a DNA sequence and returns positions and translations. |
| [`gc_content`](#gc_content) | `sequence` | Calculates GC percentage and base composition of a DNA or RNA sequence. |
| [`get_marker_genes`](#get_marker_genes) | `expression` | Returns filtered marker genes with statistics for a group from prior rank_genes_groups results. |
| [`get_string_enrichment`](#get_string_enrichment) | `genome_annotation` | Computes STRING functional enrichment (GO, KEGG and others) for a list of genes or proteins. |
| [`methylation_summary_by_chromosome`](#methylation_summary_by_chromosome) | `regulation_epigenome` | Summarizes methylation levels (mean, median, count) per chromosome for a loaded methylation table. |
| [`parse_bedgraph`](#parse_bedgraph) | `regulation_epigenome` | Parses a bedGraph/BED file of genomic intervals with scores and summarizes chromosomes and score distribution. |
| [`parse_fasta`](#parse_fasta) | `sequence` | Parses a FASTA file and returns sequence count, IDs, lengths and GC content. |
| [`parse_methylation_table`](#parse_methylation_table) | `regulation_epigenome` | Loads a tabular methylation calls file (e.g. Bismark, nanome) with chromosome, position, methylation and coverage columns. |
| [`parse_vcf`](#parse_vcf) | `variant` | Parses a VCF/VCF.gz file, optionally within a region, and returns variant counts, types and sample records. |
| [`plot_dotplot`](#plot_dotplot) | `expression` | Generates a dot plot of expression for a list of genes across cell groups. |
| [`plot_heatmap`](#plot_heatmap) | `expression` | Generates a heatmap of expression for a list of genes across cell groups. |
| [`plot_marker_genes`](#plot_marker_genes) | `expression` | Plots the top-ranked marker genes per group as a matrix plot from rank_genes_groups results. |
| [`plot_methylation_distribution`](#plot_methylation_distribution) | `regulation_epigenome` | Plots a histogram of methylation beta values, optionally split by group. |
| [`plot_volcano`](#plot_volcano) | `expression` | Generates a volcano plot of differential expression results for a group with fold-change and p-value thresholds. |
| [`rank_genes_groups`](#rank_genes_groups) | `expression` | Runs differential expression between groups (wilcoxon, t-test, logreg) and returns top marker genes per group. |
| [`reverse_complement`](#reverse_complement) | `sequence` | Returns the reverse complement of a DNA sequence with related info. |
| [`run_enrichment_analysis`](#run_enrichment_analysis) | `genome_annotation` | Runs Enrichr over-representation analysis on a gene list against GO, KEGG, Reactome, MSigDB or WikiPathways libraries. |
| [`run_gsea_preranked`](#run_gsea_preranked) | `genome_annotation` | Runs preranked gene set enrichment analysis using genes ranked from differential expression results. |
| [`search_clinvar`](#search_clinvar) | `variant` | Searches ClinVar for variants in a gene filtered by clinical significance and returns type and location. |
| [`search_ensembl_genes`](#search_ensembl_genes) | `genome_annotation` | Looks up a gene symbol in Ensembl and returns its ID, location, biotype and transcripts. |
| [`search_ncbi_gene`](#search_ncbi_gene) | `genome_annotation` | Searches NCBI Gene by symbol, name or keyword and returns IDs with summaries, chromosome and aliases. |
| [`search_ncbi_nucleotide`](#search_ncbi_nucleotide) | `sequence` | Searches the NCBI Nucleotide database and returns matching accessions and titles. |
| [`spatial_autocorrelation`](#spatial_autocorrelation) | `expression` | Computes Moran's I or Geary's C for all genes to identify spatially variable genes. |
| [`translate_sequence`](#translate_sequence) | `sequence` | Translates a DNA sequence to protein using a chosen NCBI codon table. |
| [`vcf_summary_stats`](#vcf_summary_stats) | `variant` | Computes VCF summary statistics: variant counts by type, chromosome distribution and Ti/Tv ratio. |

## `bam_stats` {#bam_stats}

Computes summary statistics (mapped/unmapped reads, references) from a BAM or SAM alignment file.

- **类别：**`sequencing_data`（Sequencing datasets, runs, experiments and read files）
- **必填参数：**`file_path`
- **关键词：**bam, sam, read alignment, mapping rate, samtools

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file_path` | string | 是 |  |

??? note "上游说明（英文原文）"

    Get summary statistics from a BAM/SAM file.

    Args:
    file_path: Path to a BAM or SAM file.

    Returns:
    Alignment statistics: mapped/unmapped reads, references, etc.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'bam_stats',
    {'file_path': '...'},
)
```

## `cpg_density` {#cpg_density}

Computes CpG observed/expected ratio in sliding windows along a DNA sequence and classifies windows as island, shore or open sea.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**cpg island, cpg obs/exp, sliding window, methylation, dna sequence

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 |  |
| `window_size` | integer | 否 | 默认 `100` |
| `step` | integer | 否 | 默认 `50` |

??? note "上游说明（英文原文）"

    Compute CpG density (obs/exp ratio) along a DNA sequence using a sliding window.

    CpG observed/expected = (count_CpG × N) / (count_C × count_G)

    Args:
    sequence: DNA sequence string.
    window_size: Window size in bp (default: 100).
    step: Step size in bp (default: 50).

    Returns:
    List of windows with CpG density and classification (island / shore / open sea).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'cpg_density',
    {'sequence': '...'},
)
```

## `export_gene_expression` {#export_gene_expression}

Extracts expression values for a list of genes across all cells of an AnnData dataset, optionally from a layer.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`dataset_name`, `genes`
- **关键词：**gene expression, expression matrix, anndata, single cell, counts

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 |  |
| `genes` | array of string | 是 |  |
| `layer` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Extract expression values for specific genes across all cells.

    Args:
    dataset_name: Name of the AnnData dataset.
    genes: List of gene names to extract.
    layer: Layer to use (None = .X, 'counts' for raw counts).

    Returns:
    Expression matrix as records (cells × genes).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'export_gene_expression',
    {'dataset_name': '...', 'genes': []},
)
```

## `fetch_ensembl_sequence` {#fetch_ensembl_sequence}

Fetches a genomic, cDNA or CDS nucleotide sequence from Ensembl for a gene or transcript ID.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`ensembl_id`
- **关键词：**ensembl, nucleotide sequence, cdna, cds, transcript

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `ensembl_id` | string | 是 |  |
| `seq_type` | string | 否 | 默认 `genomic` |

??? note "上游说明（英文原文）"

    Fetch a nucleotide sequence from Ensembl by ID.

    Args:
    ensembl_id: Ensembl gene/transcript ID (e.g. 'ENSG00000012048',
    'ENST00000357654').
    seq_type: 'genomic', 'cdna', or 'cds' (default: 'genomic').

    Returns:
    Sequence with metadata.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'fetch_ensembl_sequence',
    {'ensembl_id': '...'},
)
```

## `fetch_gene_ontology` {#fetch_gene_ontology}

Fetches GO biological process, molecular function and cellular component annotations for a gene symbol via UniProt.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`gene_symbol`
- **关键词：**gene ontology, go annotation, biological process, molecular function, cellular component

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_symbol` | string | 是 |  |
| `organism` | string | 否 | 默认 `human` |
| `max_results` | integer | 否 | 默认 `50` |

??? note "上游说明（英文原文）"

    Fetch Gene Ontology (GO) annotations for a gene via UniProt.

    Returns Biological Process, Molecular Function, and Cellular Component terms.

    Args:
    gene_symbol: Gene symbol (e.g. 'BRCA1').
    organism: Organism name (default: 'human').
    max_results: Maximum GO terms to return (default: 50).

    Returns:
    GO annotations grouped by category (BP, MF, CC).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'fetch_gene_ontology',
    {'gene_symbol': '...'},
)
```

## `fetch_ncbi_gene_details` {#fetch_ncbi_gene_details}

Fetches NCBI Gene details (name, summary, genomic location, references) for an NCBI Gene ID.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`gene_id`
- **关键词：**ncbi gene, entrez, gene summary, genomic location

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_id` | string | 是 |  |

??? note "上游说明（英文原文）"

    Fetch detailed information about a gene from NCBI Gene by its ID.

    Args:
    gene_id: NCBI Gene ID (e.g. '672' for BRCA1).

    Returns:
    Full gene details: name, summary, genomic location, references.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'fetch_ncbi_gene_details',
    {'gene_id': '...'},
)
```

## `fetch_ncbi_sequence` {#fetch_ncbi_sequence}

Downloads a nucleotide sequence from NCBI by accession in FASTA or GenBank format, with optional truncation.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`accession`
- **关键词：**ncbi nucleotide, accession, fasta, genbank, refseq

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `accession` | string | 是 |  |
| `seq_type` | string | 否 | 默认 `fasta` |
| `max_length` | integer | 否 | 默认 `50000` |

??? note "上游说明（英文原文）"

    Download a nucleotide sequence from NCBI by accession number.

    Args:
    accession: NCBI accession (e.g. 'NM_007294.4', 'NC_000017.11').
    seq_type: 'fasta' or 'gb' for GenBank format (default: 'fasta').
    max_length: Truncate sequences longer than this (default: 50000).

    Returns:
    Sequence in requested format with metadata.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'fetch_ncbi_sequence',
    {'accession': '...'},
)
```

## `find_differentially_methylated_regions` {#find_differentially_methylated_regions}

Identifies differentially methylated regions between two sample groups using a sliding genomic window and a minimum methylation difference.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`dataset_name`, `group_col`, `group_a`, `group_b`
- **关键词：**dmr, dna methylation, differential methylation, beta value, sliding window

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 |  |
| `group_col` | string | 是 |  |
| `group_a` | string | 是 |  |
| `group_b` | string | 是 |  |
| `chrom_col` | string | 否 | 默认 `chrom` |
| `pos_col` | string | 否 | 默认 `pos` |
| `meth_col` | string | 否 | 默认 `methylation` |
| `min_diff` | number | 否 | 默认 `0.2` |
| `window_size` | integer | 否 | 默认 `1000` |

??? note "上游说明（英文原文）"

    Identify differentially methylated regions (DMRs) between two groups.

    Uses a simple sliding-window approach: for each genomic window, compute
    the mean methylation difference between groups and flag regions where
    \|Δβ\| ≥ *min_diff*.

    Args:
    dataset_name: Name of the loaded methylation DataFrame.
    group_col: Column identifying sample groups.
    group_a: First group label.
    group_b: Second group label.
    chrom_col: Chromosome column.
    pos_col: Position column.
    meth_col: Methylation column.
    min_diff: Minimum absolute methylation difference (default: 0.2).
    window_size: Window size in bp (default: 1000).

    Returns:
    List of DMRs with coordinates and methylation differences.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'find_differentially_methylated_regions',
    {'dataset_name': '...',
     'group_col': '...',
     'group_a': '...',
     'group_b': '...'},
)
```

## `find_orfs` {#find_orfs}

Finds all open reading frames above a minimum length in a DNA sequence and returns positions and translations.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**orf, open reading frame, coding region, translation, dna sequence

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 |  |
| `min_length` | integer | 否 | 默认 `100` |

??? note "上游说明（英文原文）"

    Find all open reading frames (ORFs) in a DNA sequence.

    Args:
    sequence: DNA sequence string.
    min_length: Minimum ORF length in nucleotides (default: 100).

    Returns:
    List of ORFs with positions, lengths, and translated sequences.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'find_orfs',
    {'sequence': '...'},
)
```

## `gc_content` {#gc_content}

Calculates GC percentage and base composition of a DNA or RNA sequence.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**gc content, base composition, nucleotide, dna sequence

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 |  |

??? note "上游说明（英文原文）"

    Calculate GC content of a DNA/RNA sequence.

    Args:
    sequence: DNA or RNA sequence string (e.g., "ATCGATCG").

    Returns:
    GC content percentage and base composition.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'gc_content',
    {'sequence': '...'},
)
```

## `get_marker_genes` {#get_marker_genes}

Returns filtered marker genes with statistics for a group from prior rank_genes_groups results.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`dataset_name`, `group`
- **关键词：**marker genes, differential expression, log fold change, single cell

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 |  |
| `group` | string | 是 |  |
| `n_genes` | integer | 否 | 默认 `20` |
| `min_logfoldchange` | number | 否 | 默认 `0.25` |
| `max_pval` | number | 否 | 默认 `0.05` |

??? note "上游说明（英文原文）"

    Get marker genes for a specific group from a previous rank_genes_groups run.

    Args:
    dataset_name: Name of the AnnData dataset.
    group: Group name to get markers for.
    n_genes: Maximum number of genes to return (default: 20).
    min_logfoldchange: Minimum log fold change filter (default: 0.25).
    max_pval: Maximum adjusted p-value filter (default: 0.05).

    Returns:
    Filtered marker gene list with statistics.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'get_marker_genes',
    {'dataset_name': '...', 'group': '...'},
)
```

## `get_string_enrichment` {#get_string_enrichment}

Computes STRING functional enrichment (GO, KEGG and others) for a list of genes or proteins.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`genes`
- **关键词：**functional enrichment, gene set, go, kegg, string

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `genes` | array of string | 是 |  |
| `organism` | integer | 否 | 默认 `9606` |

??? note "上游说明（英文原文）"

    Get functional enrichment for a set of genes/proteins via STRING.

    Returns enriched GO terms, KEGG pathways, and other annotations.

    Args:
    genes: List of gene/protein symbols.
    organism: NCBI taxonomy ID (default: 9606 = human).

    Returns:
    Enriched functional terms grouped by category.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'get_string_enrichment',
    {'genes': []},
)
```

## `methylation_summary_by_chromosome` {#methylation_summary_by_chromosome}

Summarizes methylation levels (mean, median, count) per chromosome for a loaded methylation table.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`dataset_name`
- **关键词：**dna methylation, chromosome, summary statistics, beta value

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 |  |
| `chrom_col` | string | 否 | 默认 `chrom` |
| `meth_col` | string | 否 | 默认 `methylation` |

??? note "上游说明（英文原文）"

    Summarise methylation levels per chromosome.

    Args:
    dataset_name: Name of the loaded methylation DataFrame.
    chrom_col: Column with chromosome names.
    meth_col: Column with methylation values.

    Returns:
    Per-chromosome mean, median, and count.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'methylation_summary_by_chromosome',
    {'dataset_name': '...'},
)
```

## `parse_bedgraph` {#parse_bedgraph}

Parses a bedGraph/BED file of genomic intervals with scores and summarizes chromosomes and score distribution.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`file_path`
- **关键词：**bedgraph, bed, methylation track, genomic intervals, signal track

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file_path` | string | 是 |  |
| `name` | string \| null | 否 | 默认 `None` |
| `max_rows` | integer | 否 | 默认 `500000` |

??? note "上游说明（英文原文）"

    Parse a bedGraph file (commonly used for methylation tracks).

    Expects 4+ columns: chrom, start, end, score.  Extra columns are kept.

    Args:
    file_path: Path to a .bedGraph or .bed file.
    name: Optional dataset name.
    max_rows: Maximum rows to load (default: 500 000).

    Returns:
    Summary with shape, chromosomes, and score distribution.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'parse_bedgraph',
    {'file_path': '...'},
)
```

## `parse_fasta` {#parse_fasta}

Parses a FASTA file and returns sequence count, IDs, lengths and GC content.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`file_path`
- **关键词：**fasta, sequence parsing, gc content, sequence length

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file_path` | string | 是 |  |
| `max_sequences` | integer | 否 | 默认 `100` |

??? note "上游说明（英文原文）"

    Parse a FASTA file and return sequence summaries.

    Args:
    file_path: Path to a FASTA file.
    max_sequences: Maximum number of sequences to return (default: 100).

    Returns:
    Sequence count, IDs, lengths, and GC content.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'parse_fasta',
    {'file_path': '...'},
)
```

## `parse_methylation_table` {#parse_methylation_table}

Loads a tabular methylation calls file (e.g. Bismark, nanome) with chromosome, position, methylation and coverage columns.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`file_path`
- **关键词：**bismark, methylation calls, cpg, coverage, dna methylation

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file_path` | string | 是 |  |
| `name` | string \| null | 否 | 默认 `None` |
| `sep` | string | 否 | 默认 `` |
| `chrom_col` | string | 否 | 默认 `chrom` |
| `pos_col` | string | 否 | 默认 `pos` |
| `meth_col` | string | 否 | 默认 `methylation` |
| `coverage_col` | string | 否 | 默认 `coverage` |
| `max_rows` | integer | 否 | 默认 `1000000` |

??? note "上游说明（英文原文）"

    Load a tabular methylation calls file (e.g. Bismark, nanome output).

    Args:
    file_path: Path to TSV / CSV.
    name: Optional dataset name.
    sep: Delimiter (default: tab).
    chrom_col: Column name for chromosome.
    pos_col: Column name for position.
    meth_col: Column name for methylation level (0-1 or 0-100).
    coverage_col: Column name for coverage / read depth.
    max_rows: Maximum rows to load.

    Returns:
    Summary statistics.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'parse_methylation_table',
    {'file_path': '...'},
)
```

## `parse_vcf` {#parse_vcf}

Parses a VCF/VCF.gz file, optionally within a region, and returns variant counts, types and sample records.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`file_path`
- **关键词：**vcf, variant calls, snp, indel, genomic region

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file_path` | string | 是 |  |
| `region` | string \| null | 否 | 默认 `None` |
| `max_variants` | integer | 否 | 默认 `500` |

??? note "上游说明（英文原文）"

    Parse a VCF file and return variant summaries.

    Args:
    file_path: Path to a VCF or VCF.gz file.
    region: Genomic region to filter (e.g., "chr1:10000-20000"). Optional.
    max_variants: Maximum number of variants to return (default: 500).

    Returns:
    Variant count, types, and sample data.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'parse_vcf',
    {'file_path': '...'},
)
```

## `plot_dotplot` {#plot_dotplot}

Generates a dot plot of expression for a list of genes across cell groups.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`dataset_name`, `genes`
- **关键词：**dot plot, gene expression, visualization, scanpy

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 |  |
| `genes` | array of string | 是 |  |
| `groupby` | string | 否 | 默认 `leiden` |
| `title` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Generate a dot plot of gene expression across groups.

    Args:
    dataset_name: Name of the AnnData dataset.
    genes: List of gene names to plot.
    groupby: .obs column defining groups (default: "leiden").
    title: Optional plot title.

    Returns:
    Base64-encoded PNG image of the dot plot.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'plot_dotplot',
    {'dataset_name': '...', 'genes': []},
)
```

## `plot_heatmap` {#plot_heatmap}

Generates a heatmap of expression for a list of genes across cell groups.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`dataset_name`, `genes`
- **关键词：**heatmap, gene expression, visualization, scanpy

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 |  |
| `genes` | array of string | 是 |  |
| `groupby` | string | 否 | 默认 `leiden` |
| `title` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Generate a heatmap of gene expression across groups.

    Args:
    dataset_name: Name of the AnnData dataset.
    genes: List of gene names to plot.
    groupby: .obs column defining groups (default: "leiden").
    title: Optional plot title.

    Returns:
    Base64-encoded PNG image of the heatmap.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'plot_heatmap',
    {'dataset_name': '...', 'genes': []},
)
```

## `plot_marker_genes` {#plot_marker_genes}

Plots the top-ranked marker genes per group as a matrix plot from rank_genes_groups results.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`dataset_name`
- **关键词：**marker genes, matrix plot, differential expression, visualization

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 |  |
| `n_genes` | integer | 否 | 默认 `5` |
| `groupby` | string | 否 | 默认 `leiden` |

??? note "上游说明（英文原文）"

    Plot top marker genes per group as a ranked genes matrix plot.

    Requires rank_genes_groups to have been run first.

    Args:
    dataset_name: Name of the AnnData dataset.
    n_genes: Number of top genes per group to display (default: 5).
    groupby: .obs column that was used for DE (default: "leiden").

    Returns:
    Base64-encoded PNG image.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'plot_marker_genes',
    {'dataset_name': '...'},
)
```

## `plot_methylation_distribution` {#plot_methylation_distribution}

Plots a histogram of methylation beta values, optionally split by group.

- **类别：**`regulation_epigenome`（Regulatory elements, TF binding, chromatin, methylation and 3D genome）
- **必填参数：**`dataset_name`
- **关键词：**dna methylation, beta value, histogram, visualization

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 |  |
| `meth_col` | string | 否 | 默认 `methylation` |
| `group_col` | string \| null | 否 | 默认 `None` |

??? note "上游说明（英文原文）"

    Plot the distribution of methylation values (beta-value histogram).

    Args:
    dataset_name: Name of the loaded methylation DataFrame.
    meth_col: Column with methylation values.
    group_col: Optional column to split by (overlay histograms).

    Returns:
    Base64-encoded PNG histogram.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'plot_methylation_distribution',
    {'dataset_name': '...'},
)
```

## `plot_volcano` {#plot_volcano}

Generates a volcano plot of differential expression results for a group with fold-change and p-value thresholds.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`dataset_name`, `group`
- **关键词：**volcano plot, differential expression, fold change, p-value

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 |  |
| `group` | string | 是 |  |
| `log2fc_threshold` | number | 否 | 默认 `1.0` |
| `pval_threshold` | number | 否 | 默认 `0.05` |
| `top_n_label` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Generate a volcano plot from differential expression results.

    Requires rank_genes_groups to have been run first.

    Args:
    dataset_name: Name of the AnnData dataset.
    group: Group from DE results to plot.
    log2fc_threshold: Log2 fold change threshold for significance lines (default: 1.0).
    pval_threshold: Adjusted p-value threshold (default: 0.05).
    top_n_label: Number of top genes to label (default: 10).

    Returns:
    Base64-encoded PNG image of the volcano plot.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'plot_volcano',
    {'dataset_name': '...', 'group': '...'},
)
```

## `rank_genes_groups` {#rank_genes_groups}

Runs differential expression between groups (wilcoxon, t-test, logreg) and returns top marker genes per group.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`dataset_name`
- **关键词：**differential expression, marker genes, wilcoxon, single cell

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 |  |
| `groupby` | string | 否 | 默认 `leiden` |
| `method` | string | 否 | 默认 `wilcoxon` |
| `n_genes` | integer | 否 | 默认 `100` |
| `reference` | string | 否 | 默认 `rest` |

??? note "上游说明（英文原文）"

    Run differential expression analysis between groups.

    Args:
    dataset_name: Name of the AnnData dataset.
    groupby: .obs column defining groups (default: "leiden").
    method: Test method — 'wilcoxon', 't-test', 't-test_overestim_var',
    'logreg' (default: 'wilcoxon').
    n_genes: Number of top genes to return per group (default: 100).
    reference: Reference group — 'rest' for one-vs-rest (default: 'rest').

    Returns:
    Top marker genes per group.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'rank_genes_groups',
    {'dataset_name': '...'},
)
```

## `reverse_complement` {#reverse_complement}

Returns the reverse complement of a DNA sequence with related info.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**reverse complement, dna strand, complement, nucleotide

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 |  |

??? note "上游说明（英文原文）"

    Get the reverse complement of a DNA sequence.

    Args:
    sequence: DNA sequence string.

    Returns:
    Reverse complement and related info.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'reverse_complement',
    {'sequence': '...'},
)
```

## `run_enrichment_analysis` {#run_enrichment_analysis}

Runs Enrichr over-representation analysis on a gene list against GO, KEGG, Reactome, MSigDB or WikiPathways libraries.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`gene_list`
- **关键词：**enrichr, over-representation, gene set, pathway enrichment, go

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_list` | array of string | 是 |  |
| `gene_set_library` | string | 否 | 默认 `GO_Biological_Process_2023` |
| `organism` | string | 否 | 默认 `human` |
| `top_n` | integer | 否 | 默认 `20` |

??? note "上游说明（英文原文）"

    Run over-representation analysis (ORA) on a gene list using Enrichr libraries.

    Args:
    gene_list: List of gene symbols to test.
    gene_set_library: Enrichr library name. Common choices:
    - GO_Biological_Process_2023
    - GO_Molecular_Function_2023
    - KEGG_2021_Human
    - Reactome_2022
    - MSigDB_Hallmark_2020
    - WikiPathway_2023_Human
    organism: 'human' or 'mouse' (default: 'human').
    top_n: Number of top results to return (default: 20).

    Returns:
    Enrichment results with terms, p-values, and gene overlaps.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'run_enrichment_analysis',
    {'gene_list': []},
)
```

## `run_gsea_preranked` {#run_gsea_preranked}

Runs preranked gene set enrichment analysis using genes ranked from differential expression results.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`dataset_name`, `group`
- **关键词：**gsea, gene set enrichment, msigdb, ranked genes

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 |  |
| `group` | string | 是 |  |
| `gene_set_library` | string | 否 | 默认 `GO_Biological_Process_2023` |
| `top_n` | integer | 否 | 默认 `20` |

??? note "上游说明（英文原文）"

    Run Gene Set Enrichment Analysis (GSEA) using pre-ranked genes from DE results.

    Requires rank_genes_groups to have been run first.

    Args:
    dataset_name: Name of the AnnData dataset (must have DE results).
    group: Group from the DE analysis to use for ranking.
    gene_set_library: Enrichr/MSigDB library name.
    top_n: Number of top results to return.

    Returns:
    GSEA results with normalized enrichment scores.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'run_gsea_preranked',
    {'dataset_name': '...', 'group': '...'},
)
```

## `search_clinvar` {#search_clinvar}

Searches ClinVar for variants in a gene filtered by clinical significance and returns type and location.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`gene`
- **关键词：**clinvar, pathogenic, clinical significance, gene variants

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene` | string | 是 |  |
| `significance` | string | 否 | 默认 `pathogenic` |
| `max_results` | integer | 否 | 默认 `20` |

??? note "上游说明（英文原文）"

    Search ClinVar for clinical variants associated with a gene.

    Args:
    gene: Gene symbol (e.g. 'BRCA1', 'TP53').
    significance: Clinical significance filter — 'pathogenic',
    'likely_pathogenic', 'benign', 'uncertain', or 'all'
    (default: 'pathogenic').
    max_results: Maximum results (default: 20).

    Returns:
    List of clinical variants with significance, variant type, and location.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'search_clinvar',
    {'gene': '...'},
)
```

## `search_ensembl_genes` {#search_ensembl_genes}

Looks up a gene symbol in Ensembl and returns its ID, location, biotype and transcripts.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`symbol`
- **关键词：**ensembl, gene symbol, ensembl id, transcripts, biotype

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `symbol` | string | 是 |  |
| `species` | string | 否 | 默认 `homo_sapiens` |

??? note "上游说明（英文原文）"

    Look up a gene in Ensembl by symbol.

    Args:
    symbol: Gene symbol (e.g. 'BRCA1', 'TP53').
    species: Ensembl species name (default: 'homo_sapiens').

    Returns:
    Gene info with Ensembl ID, location, biotype, and transcripts.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'search_ensembl_genes',
    {'symbol': '...'},
)
```

## `search_ncbi_gene` {#search_ncbi_gene}

Searches NCBI Gene by symbol, name or keyword and returns IDs with summaries, chromosome and aliases.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`query`
- **关键词：**ncbi gene, entrez, gene search, aliases

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 |  |
| `organism` | string | 否 | 默认 `human` |
| `max_results` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Search NCBI Gene database by gene symbol, name, or keyword.

    Args:
    query: Gene symbol (e.g. 'BRCA1'), name, or free-text query.
    organism: Organism filter (default: 'human'). Use common or
    scientific names.
    max_results: Maximum number of results (default: 10).

    Returns:
    Gene IDs with summaries (name, description, chromosome, aliases).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'search_ncbi_gene',
    {'query': '...'},
)
```

## `search_ncbi_nucleotide` {#search_ncbi_nucleotide}

Searches the NCBI Nucleotide database and returns matching accessions and titles.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`query`
- **关键词：**ncbi nucleotide, accession, mrna, sequence search

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `query` | string | 是 |  |
| `max_results` | integer | 否 | 默认 `10` |

??? note "上游说明（英文原文）"

    Search NCBI Nucleotide database.

    Args:
    query: Search query (e.g. 'BRCA1 mRNA human', 'NC_000017.11').
    max_results: Maximum results to return (default: 10).

    Returns:
    List of matching accessions with titles.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'search_ncbi_nucleotide',
    {'query': '...'},
)
```

## `spatial_autocorrelation` {#spatial_autocorrelation}

Computes Moran's I or Geary's C for all genes to identify spatially variable genes.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`dataset_name`
- **关键词：**spatially variable genes, moran's i, geary's c, spatial transcriptomics

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `dataset_name` | string | 是 |  |
| `mode` | string | 否 | 默认 `moran` |
| `n_jobs` | integer | 否 | 默认 `1` |
| `n_perms` | integer | 否 | 默认 `100` |

??? note "上游说明（英文原文）"

    Compute spatial autocorrelation statistics for all genes.

    Identifies spatially variable genes using Moran's I or Geary's C.

    Args:
    dataset_name: Name of the AnnData dataset.
    mode: 'moran' for Moran's I or 'geary' for Geary's C (default: 'moran').
    n_jobs: Number of parallel jobs (default: 1).
    n_perms: Number of permutations for p-value (default: 100).

    Returns:
    Top spatially variable genes with statistics.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'spatial_autocorrelation',
    {'dataset_name': '...'},
)
```

## `translate_sequence` {#translate_sequence}

Translates a DNA sequence to protein using a chosen NCBI codon table.

- **类别：**`sequence`（DNA/RNA sequence retrieval, analysis, alignment, simulation or manipulation）
- **必填参数：**`sequence`
- **关键词：**translation, codon table, protein sequence, dna sequence

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 |  |
| `table` | integer | 否 | 默认 `1` |

??? note "上游说明（英文原文）"

    Translate a DNA sequence to protein.

    Args:
    sequence: DNA sequence (length should be divisible by 3).
    table: NCBI codon table number (default: 1 = Standard).

    Returns:
    Protein sequence and properties.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'translate_sequence',
    {'sequence': '...'},
)
```

## `vcf_summary_stats` {#vcf_summary_stats}

Computes VCF summary statistics: variant counts by type, chromosome distribution and Ti/Tv ratio.

- **类别：**`variant`（Genetic variants, genotypes, allele frequencies, variant effect and interpretation, polygenic scores）
- **必填参数：**`file_path`
- **关键词：**vcf, ti/tv, variant counts, chromosome distribution

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `file_path` | string | 是 |  |

??? note "上游说明（英文原文）"

    Calculate summary statistics from a VCF file.

    Args:
    file_path: Path to a VCF or VCF.gz file.

    Returns:
    Summary: variant counts by type, chromosome distribution, Ti/Tv ratio.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'bioinformatics-mcp',
    'vcf_summary_stats',
    {'file_path': '...'},
)
```
