# ToolUniverse: regulation_epigenome (`tooluniverse`)

DNAKit's MCP tool index keeps **845** DNA/gene tools of this server (2927 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

This page lists the 54 tools in category `regulation_epigenome`; other categories: [`sequence`](tooluniverse-sequence.md) (44), [`genome_annotation`](tooluniverse-genome_annotation.md) (253), [`variant`](tooluniverse-variant.md) (155), [`association`](tooluniverse-association.md) (138), [`expression`](tooluniverse-expression.md) (56), [`sequencing_data`](tooluniverse-sequencing_data.md) (91), [`prediction_model`](tooluniverse-prediction_model.md) (17), [`molecular_tools`](tooluniverse-molecular_tools.md) (23), [`breeding_population`](tooluniverse-breeding_population.md) (14).

| Item | Value |
| --- | --- |
| Source | [mims-harvard/ToolUniverse @ `9959ebe`](https://github.com/mims-harvard/ToolUniverse/tree/9959ebee7638773314362f1746980a15f5b17e5b) |
| License | Apache-2.0 |
| Command | `tooluniverse` |
| Transports | stdio |
| Installation | Install the pinned ToolUniverse distribution and required provider dependencies separately. Use the installed tooluniverse console script; unlike the README uvx --refresh example, this bridge does not refresh or install anything. |
| Effects | Starts the compact ToolUniverse MCP catalog. Selected tools may contact services, load models, execute programs or write files according to their provider definitions. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`Activity_infer_ulm`](#Activity_infer_ulm) | `regulation_epigenome` | Infers TF or pathway activities from per-gene statistics and a regulator-target network using decoupleR ULM. |
| [`ChIPAtlas_enrichment_analysis`](#ChIPAtlas_enrichment_analysis) | `regulation_epigenome` | Runs ChIP-Atlas enrichment of transcription factors and histone marks for BED regions, motifs or genes. |
| [`ChIPAtlas_get_colocalization`](#ChIPAtlas_get_colocalization) | `regulation_epigenome` | Returns proteins whose ChIP-seq peaks colocalize with a given antigen in a cell type class. |
| [`ChIPAtlas_get_peak_data`](#ChIPAtlas_get_peak_data) | `regulation_epigenome` | Returns download URLs for ChIP-Atlas peak calls and coverage in BigWig, BED or BigBed format. |
| [`ChIPAtlas_get_target_genes`](#ChIPAtlas_get_target_genes) | `regulation_epigenome` | Returns a TF's ranked potential target genes based on TSS-proximal ChIP-seq binding. |
| [`ENCODE_get_chromatin_state`](#ENCODE_get_chromatin_state) | `regulation_epigenome` | Searches ENCODE ChromHMM chromatin state annotations by biosample. |
| [`ENCODE_search_annotations`](#ENCODE_search_annotations) | `regulation_epigenome` | Searches ENCODE annotations such as cCREs, chromatin states and imputed signals. |
| [`ENCODE_search_chromatin_accessibility`](#ENCODE_search_chromatin_accessibility) | `regulation_epigenome` | Searches ENCODE ATAC-seq and DNase-seq experiments. |
| [`ENCODE_search_hic_experiments`](#ENCODE_search_hic_experiments) | `regulation_epigenome` | Searches ENCODE Hi-C experiments. |
| [`ENCODE_search_histone_experiments`](#ENCODE_search_histone_experiments) | `regulation_epigenome` | Searches ENCODE histone ChIP-seq experiments by mark and biosample. |
| [`ENCODE_search_methylation_experiments`](#ENCODE_search_methylation_experiments) | `regulation_epigenome` | Searches ENCODE WGBS and RRBS methylation experiments. |
| [`ENCORI_get_ceRNA_network`](#ENCORI_get_ceRNA_network) | `regulation_epigenome` | Builds ENCORI ceRNA networks of genes and lncRNAs sharing miRNA binding with a query gene. |
| [`ENCORI_get_degradome_cleavage`](#ENCORI_get_degradome_cleavage) | `regulation_epigenome` | Retrieves ENCORI degradome-seq validated miRNA cleavage sites by gene or miRNA. |
| [`ENCORI_get_miRNA_targets`](#ENCORI_get_miRNA_targets) | `regulation_epigenome` | Looks up ENCORI CLIP-supported and predicted miRNA-target interactions by miRNA or gene. |
| [`ENCORI_get_RBP_disease`](#ENCORI_get_RBP_disease) | `regulation_epigenome` | Retrieves ENCORI RBP binding sites on genes that overlap COSMIC somatic mutations by tissue or disease. |
| [`ENCORI_get_RBP_targets`](#ENCORI_get_RBP_targets) | `regulation_epigenome` | Looks up ENCORI CLIP-seq RBP binding sites on transcripts, by RBP or target gene. |
| [`ENCORI_get_RNA_RNA_interactions`](#ENCORI_get_RNA_RNA_interactions) | `regulation_epigenome` | Retrieves ENCORI RNA-RNA interaction partners of an RNA from PARIS/LIGR-seq/SPLASH data. |
| [`ENCORI_scan_RBP_motifs`](#ENCORI_scan_RBP_motifs) | `regulation_epigenome` | Scans ENCORI RBP binding motif enrichment by motif sequence or RBP. |
| [`ensembl_get_regulatory_features`](#ensembl_get_regulatory_features) | `regulation_epigenome` | Gets Ensembl regulatory features such as promoters, enhancers and TF binding sites in a region. |
| [`EnsemblReg_get_binding_matrix`](#EnsemblReg_get_binding_matrix) | `regulation_epigenome` | Gets an Ensembl TF binding matrix by stable ID with frequencies and consensus. |
| [`EnsemblReg_get_motif_features`](#EnsemblReg_get_motif_features) | `regulation_epigenome` | Gets TF binding motif features in a genomic region from Ensembl. |
| [`EnsemblReg_get_regulatory_elements`](#EnsemblReg_get_regulatory_elements) | `regulation_epigenome` | Gets Ensembl Regulatory Build features (enhancers, promoters, CTCF sites, open chromatin) for a region. |
| [`EWASCatalog_search_by_cpg`](#EWASCatalog_search_by_cpg) | `regulation_epigenome` | Gets EWAS Catalog associations for a CpG methylation site. |
| [`EWASCatalog_search_by_gene`](#EWASCatalog_search_by_gene) | `regulation_epigenome` | Gets EWAS Catalog associations across CpG sites mapped to a gene. |
| [`FourDN_get_file_metadata`](#FourDN_get_file_metadata) | `regulation_epigenome` | Gets 4DN file metadata for Hi-C contact matrices, TAD calls, A/B compartments and chromatin loops. |
| [`FourDN_search_data`](#FourDN_search_data) | `regulation_epigenome` | Searches the 4DN Data Portal for Hi-C chromatin conformation files and experiments by assay, cell type or keyword. |
| [`GEO_search_atacseq_datasets`](#GEO_search_atacseq_datasets) | `regulation_epigenome` | Searches NCBI GEO for ATAC-seq datasets. |
| [`GEO_search_chipseq_datasets`](#GEO_search_chipseq_datasets) | `regulation_epigenome` | Searches NCBI GEO for ChIP-seq datasets. |
| [`GEO_search_methylation_datasets`](#GEO_search_methylation_datasets) | `regulation_epigenome` | Searches NCBI GEO for DNA methylation array datasets (450K, EPIC). |
| [`HOCOMOCO_get_motif`](#HOCOMOCO_get_motif) | `regulation_epigenome` | Gets a HOCOMOCO TF binding motif with gene names, TF class, quality, consensus and optional PWM. |
| [`HOCOMOCO_search_motifs`](#HOCOMOCO_search_motifs) | `regulation_epigenome` | Searches HOCOMOCO v14 TF binding motifs by gene/protein name returning IDs, quality and consensus. |
| [`jaspar_get_matrix`](#jaspar_get_matrix) | `regulation_epigenome` | Returns metadata and the PFM/PWM for a single JASPAR matrix by matrix ID. |
| [`jaspar_get_matrix_versions`](#jaspar_get_matrix_versions) | `regulation_epigenome` | Lists all versions of a JASPAR TF binding matrix for a base ID. |
| [`JASPAR_get_transcription_factors`](#JASPAR_get_transcription_factors) | `regulation_epigenome` | Lists all JASPAR transcription factor binding profile matrices with pagination and optional collection filter. |
| [`jaspar_search_matrices`](#jaspar_search_matrices) | `regulation_epigenome` | Searches and filters JASPAR TF binding matrices by text, TF name, collection, taxon or species. |
| [`MEME_discover_motifs`](#MEME_discover_motifs) | `regulation_epigenome` | Runs MEME de novo motif discovery on a set of DNA sequences. |
| [`MEME_fimo_scan`](#MEME_fimo_scan) | `regulation_epigenome` | Scans DNA sequences for occurrences of a MEME-format motif using FIMO and returns significant sites. |
| [`MEME_tomtom_compare`](#MEME_tomtom_compare) | `regulation_epigenome` | Compares a query motif against known TF motif databases with TOMTOM and returns similar motifs. |
| [`OmniPath_get_dorothea_regulon`](#OmniPath_get_dorothea_regulon) | `regulation_epigenome` | Returns the DoRothEA regulon of a transcription factor with target genes, mode of regulation and confidence. |
| [`OmniPath_get_tf_target_interactions`](#OmniPath_get_tf_target_interactions) | `regulation_epigenome` | Returns TF-target gene interactions from DoRothEA and CollecTRI with mode of regulation and confidence. |
| [`RegulomeDB_query_variant`](#RegulomeDB_query_variant) | `regulation_epigenome` | Returns RegulomeDB regulatory evidence and score for a variant by rsID. |
| [`ReMap_get_peaks_in_region`](#ReMap_get_peaks_in_region) | `regulation_epigenome` | Returns ReMap ChIP-seq transcriptional regulator peaks overlapping a genomic interval. |
| [`ReMap_get_transcription_factor_binding`](#ReMap_get_transcription_factor_binding) | `regulation_epigenome` | Returns ReMap transcription factor binding sites for a gene and cell type. |
| [`ReMap_list_datasets_for_target`](#ReMap_list_datasets_for_target) | `regulation_epigenome` | Lists ReMap ChIP-seq datasets for a transcription factor with optional peak counts. |
| [`SCREEN_get_cCRE_by_accession`](#SCREEN_get_cCRE_by_accession) | `regulation_epigenome` | Looks up ENCODE SCREEN candidate cis-regulatory elements by accession and returns their classification (promoter-like, enhancer-like, CTCF) and genomic locatio… |
| [`SCREEN_get_regulatory_elements`](#SCREEN_get_regulatory_elements) | `regulation_epigenome` | Retrieves ENCODE SCREEN candidate cis-regulatory elements (enhancer, promoter or insulator types) associated with a gene symbol. |
| [`SCREEN_search_cCREs_by_region`](#SCREEN_search_cCREs_by_region) | `regulation_epigenome` | Queries the ENCODE SCREEN registry for candidate cis-regulatory elements overlapping a chromosome region, optionally filtered by element type. |
| [`SGD_get_regulation`](#SGD_get_regulation) | `regulation_epigenome` | Gets the curated transcriptional regulation network for a yeast gene from SGD (regulator-target pairs with evidence and direction). |
| [`UCSC_get_cpg_islands`](#UCSC_get_cpg_islands) | `regulation_epigenome` | Gets UCSC CpG island annotations for a genomic region. |
| [`UCSC_get_encode_cCREs`](#UCSC_get_encode_cCREs) | `regulation_epigenome` | Gets ENCODE4 candidate cis-regulatory elements for a region from UCSC. |
| [`UCSC_get_tf_binding_clusters`](#UCSC_get_tf_binding_clusters) | `regulation_epigenome` | Gets ENCODE3 TF ChIP-seq clusters for a region from UCSC. |
| [`UniBind_get_dataset`](#UniBind_get_dataset) | `regulation_epigenome` | Retrieves UniBind direct transcription factor-DNA binding-site details for one dataset, including TF, cell type and motif information. |
| [`UniBind_list_tfs`](#UniBind_list_tfs) | `regulation_epigenome` | Lists transcription factors profiled in UniBind, optionally filtered by TF gene symbol substring. |
| [`UniBind_search_datasets`](#UniBind_search_datasets) | `regulation_epigenome` | Searches UniBind TF-DNA binding-site datasets (one ChIP-seq experiment per dataset) by TF symbol, cell type or species. |

## `Activity_infer_ulm` {#Activity_infer_ulm}

Infers TF or pathway activities from per-gene statistics and a regulator-target network using decoupleR ULM.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `network`
- **Keywords:** decoupler, tf activity, regulon, ulm, differential expression
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_stats` | object \| null | no | \{gene: statistic\} mapping (e.g. DE t-statistic per gene). Alternative to genes+stats. |
| `genes` | array \| null | no | Gene identifiers (use with `stats`, same order). |
| `stats` | array \| null | no | Per-gene statistic (same order as `genes`). |
| `network` | array of object | yes | Regulatory edges: \[\{source, target, weight\}\] (weight defaults to 1; 'mor' accepted as an alias). |
| `min_targets` | integer \| null | no | Minimum targets present in the gene list for a source to be tested (default 5). |

??? note "Upstream description"

    Infer transcription-factor / pathway activities from a per-gene statistic (e.g. a DESeq2/edgeR t-stat or logFC) and a regulatory network (source=TF/pathway, target=gene, weight=mode of regulation), using decoupleR's univariate-linear-model (ULM) method. Returns per-source activity (slope t-value; +=active) and p-value. Pure-compute; pairs with the DoRothEA/PROGENy networks from OmniPath_get_dorothea_regulon.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Activity_infer_ulm', 'arguments': {'network': []}},
)
```

## `ChIPAtlas_enrichment_analysis` {#ChIPAtlas_enrichment_analysis}

Runs ChIP-Atlas enrichment of transcription factors and histone marks for BED regions, motifs or genes.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** chip-seq, tf enrichment, histone, bed regions
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | One of `enrichment_analysis` Default `enrichment_analysis` |
| `bed_data` | string | no | **Option 1**: BED format genomic regions (tab-separated: chr, start, end). For finding proteins bound to specific genomic regions. Example: 'chr1\\t1000\\t2000\\nchr2\\t3000\\t4000'. |
| `motif` | string | no | **Option 2**: DNA sequence motif in IUPAC notation. Use: A/T/G/C (bases), W=A\|T, S=G\|C, M=A\|C, K=G\|T, R=A\|G, Y=C\|T, B=C\|G\|T, D=A\|G\|T, H=A\|C\|T, V=A\|C\|G, N=any. For finding proteins bound to specific DNA sequences. Example: 'CANNTG' (E-box motif). |
| `gene_list` | array \| string | no | **Option 3**: Gene symbols (HGNC for human, MGI for mouse, RGD for rat, FlyBase, WormBase, SGD for yeast). Provide as array or single gene. For finding transcription factors regulating genes. Example: \['TP53', 'MDM2', 'CDKN1A'\]. |
| `genome` | string | no | Genome assembly One of `hg38`, `hg19`, `mm10`, `mm9`, `rn6`, `dm6`, `dm3`, `ce11`, `ce10`, `sacCer3` Default `hg38` |
| `antigen_class` | string | no | Filter by antigen class (e.g., 'TFs and others', 'Histone', 'RNA polymerase') |
| `cell_type_class` | string | no | Filter by cell type class (e.g., 'Blood', 'Liver', 'Brain') |
| `threshold` | string | no | Peak calling stringency (MACS2 Q-value). Options: '05'=1e-5 (permissive, more peaks, broader features), '10'=1e-10 (moderate, balanced), '20'=1e-20 (strict, high confidence only, narrow peaks). Default '05' suitable for most analyses. Higher values = fewer but more confident peaks. One of `05`, `10`, `20` Default `05` |
| `distance` | string | no | Distance from Transcription Start Site (TSS) in base pairs for gene-TF association. Defines promoter region. Default 5000 (±5kb, captures typical promoters). Use 1000-2000 for narrow promoters, 10000+ for enhancer regions. Default `5000` |

??? note "Upstream description"

    Perform enrichment analysis to identify transcription factors and histone modifications enriched in your data. **Required: Provide ONE input type** - (1) BED genomic regions, (2) DNA sequence motif (IUPAC notation: A/T/G/C/W/S/M/K/R/Y/B/D/H/V/N), or (3) gene symbol list. Compares your input against 433,000+ ChIP-seq/ATAC-seq/Bisulfite-seq experiments to identify significant enrichment. Returns ranked list of proteins bound to your regions/motif or regulating your genes. Note: Returns submission URL (web form-based analysis). Use for: identifying regulators of genomic regions, finding proteins bound to sequence motifs, discovering transcription factors regulating gene sets.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ChIPAtlas_enrichment_analysis', 'arguments': {}},
)
```

## `ChIPAtlas_get_colocalization` {#ChIPAtlas_get_colocalization}

Returns proteins whose ChIP-seq peaks colocalize with a given antigen in a cell type class.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `antigen`, `cell_type_class`
- **Keywords:** colocalization, co-binding, transcription factor, chip-seq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | One of `get_colocalization` Default `get_colocalization` |
| `antigen` | string | yes | Antigen / transcription factor symbol (e.g. 'CTCF', 'AFF4'). |
| `cell_type_class` | string | yes | Tissue / cell-type class for which the colocalization matrix was computed (e.g. 'Blood', 'Breast', 'Liver', 'Digestive tract'). Spaces are allowed. |
| `genome` | string | no | Genome assembly. One of `hg38`, `hg19`, `mm10`, `mm9`, `rn6`, `dm6`, `dm3`, `ce11`, `ce10`, `sacCer3` Default `hg38` |
| `limit` | integer | no | Maximum number of ranked partner proteins to return. Default `50` |

??? note "Upstream description"

    ChIP-Atlas Colocalization: for a given antigen/transcription factor in a specific tissue class, return the other proteins whose ChIP-seq peaks co-occur (co-bind) genome-wide, ranked by overlap score. Reveals protein complexes and co-regulators (e.g. CTCF in Blood co-binds the cohesin subunits RAD21, SMC1A, SMC3, STAG1 and ZNF143). Requires both an antigen and a cell-type/tissue class because the matrix is precomputed per antigen-per-tissue. To discover which tissue classes exist for an antigen, consult the index at https://chip-atlas.org/data/colo_analysis.json?genome=&lt;genome&gt;. Use for: identifying co-binding partners, inferring TF complexes, finding co-regulators of a factor in a tissue.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ChIPAtlas_get_colocalization',
     'arguments': {'antigen': '...', 'cell_type_class': '...'}},
)
```

## `ChIPAtlas_get_peak_data` {#ChIPAtlas_get_peak_data}

Returns download URLs for ChIP-Atlas peak calls and coverage in BigWig, BED or BigBed format.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `experiment_id`
- **Keywords:** peak calls, bigwig, bed, macs2, chip-seq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | One of `get_peak_data` Default `get_peak_data` |
| `experiment_id` | string | yes | Experiment ID (SRX/ERX/DRX format, required) |
| `genome` | string | no | Genome assembly One of `hg38`, `hg19`, `mm10`, `mm9`, `rn6`, `dm6`, `dm3`, `ce11`, `ce10`, `sacCer3` Default `hg38` |
| `format` | string | no | Output format One of `bigwig`, `bed`, `bigbed` Default `bigwig` |
| `threshold` | string | no | Q-value threshold for BED/BigBed peak files. '05'=1e-5 (more peaks), '10'=1e-10 (moderate), '20'=1e-20 (high confidence only). Only applies to BED/BigBed formats. Default '05'. One of `05`, `10`, `20` Default `05` |

??? note "Upstream description"

    Get download URLs for ChIP-Atlas peak-call data in BigWig, BED, or BigBed format. BigWig contains coverage scores in RPM (Reads Per Million). BED/BigBed contain peak regions with MACS2 (Model-based Analysis of ChIP-Seq) peak caller scores. Use for: downloading raw data, integrating with genome browsers (UCSC, IGV), custom analysis, visualization.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ChIPAtlas_get_peak_data',
     'arguments': {'experiment_id': '...'}},
)
```

## `ChIPAtlas_get_target_genes` {#ChIPAtlas_get_target_genes}

Returns a TF's ranked potential target genes based on TSS-proximal ChIP-seq binding.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `antigen`
- **Keywords:** tf target genes, tss binding, transcription factor, chip-seq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | One of `get_target_genes` Default `get_target_genes` |
| `antigen` | string | yes | Transcription factor / antigen symbol (e.g. 'GATA1', 'CTCF'). |
| `genome` | string | no | Genome assembly. One of `hg38`, `hg19`, `mm10`, `mm9`, `rn6`, `dm6`, `dm3`, `ce11`, `ce10`, `sacCer3` Default `hg38` |
| `distance` | string | no | TSS-distance window in kb defining 'bound near a gene'. '1' (±1kb, strict promoter), '5' (±5kb, default), '10' (±10kb, includes nearby enhancers). One of `1`, `5`, `10` Default `5` |
| `limit` | integer | no | Maximum number of top-scoring target genes to return. Default `100` |

??? note "Upstream description"

    ChIP-Atlas Target Genes: for a transcription factor, return its ranked list of potential target genes - genes whose TSS-proximal region is bound by the TF - scored as the average binding strength across all ChIP-seq experiments for that TF. Choose the TSS-distance window with the 'distance' parameter (1, 5, or 10 kb). Use for: predicting which genes a TF regulates, building TF-target gene sets, prioritizing direct targets. To check which TFs have target-gene data for a genome, see https://chip-atlas.org/data/target_genes_analysis.json. Note: a few antigens (including CTCF on hg38) are served with a header but no scored rows; in that case the tool returns gene_count 0 with an explanatory note - pick another factor (e.g. 'GATA1').

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ChIPAtlas_get_target_genes',
     'arguments': {'antigen': '...'}},
)
```

## `ENCODE_get_chromatin_state` {#ENCODE_get_chromatin_state}

Searches ENCODE ChromHMM chromatin state annotations by biosample.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** chromhmm, chromatin state, encode
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `biosample_term_name` | string \| null | no | Biosample name (e.g., 'K562', 'HepG2', 'GM12878', 'liver'). Leave empty for all. |
| `organism` | string | no | Organism scientific name. Default `Homo sapiens` |
| `limit` | integer | no | Maximum number of results. Default `25` |

??? note "Upstream description"

    Search ENCODE chromatin state annotations (ChromHMM segmentations) for specific biosamples or tissues. Chromatin states are genome-wide annotations that classify each genomic position into a functional state (e.g., active TSS, flanking active TSS, strong transcription, weak transcription, genic enhancers, active enhancers, bivalent enhancer, bivalent/poised TSS, flanking bivalent TSS, repressed polycomb, weak repressed polycomb, quiescent, heterochromatin) based on combinatorial patterns of histone marks. Use this to understand the chromatin landscape of specific cell types.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_get_chromatin_state', 'arguments': {}},
)
```

## `ENCODE_search_annotations` {#ENCODE_search_annotations}

Searches ENCODE annotations such as cCREs, chromatin states and imputed signals.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** ccre, enhancer, promoter, encode
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `annotation_type` | string | no | Annotation type filter. Options: 'candidate Cis-Regulatory Elements' (cCREs), 'chromatin state' (ChromHMM), 'imputed signal' (predicted epigenetic marks). Default: cCREs. Default `candidate Cis-Regulatory Elements` |
| `biosample_term_name` | string \| null | no | Biosample filter (e.g., 'K562', 'HepG2'). Leave empty for all. |
| `organism` | string | no | Organism scientific name. Default `Homo sapiens` |
| `assembly` | string | no | Genome assembly (e.g., 'GRCh38', 'hg19', 'mm10'). Default `GRCh38` |
| `limit` | integer | no | Maximum number of results. Default `25` |

??? note "Upstream description"

    Search ENCODE annotations including candidate cis-Regulatory Elements (cCREs), chromatin states, and imputed signals. cCREs are classified as promoter-like (PLS), proximal enhancer-like (pELS), distal enhancer-like (dELS), DNase-H3K4me3 (DNase-H3K4me3), and CTCF-only. Chromatin state annotations use ChromHMM models to segment the genome into functional states. Use this to find regulatory element annotations for specific biosamples or genome assemblies.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_search_annotations', 'arguments': {}},
)
```

## `ENCODE_search_chromatin_accessibility` {#ENCODE_search_chromatin_accessibility}

Searches ENCODE ATAC-seq and DNase-seq experiments.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** atac-seq, dnase-seq, open chromatin, encode
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `assay_type` | string | no | Chromatin accessibility assay: 'ATAC-seq' (Assay for Transposase-Accessible Chromatin, newer/preferred) or 'DNase-seq' (DNase I hypersensitive sites, larger dataset). Default `ATAC-seq` |
| `biosample_term_name` | string \| null | no | Biosample name filter (e.g., 'K562', 'GM12878', 'liver'). Leave empty for all. |
| `organism` | string | no | Organism scientific name. Default `Homo sapiens` |
| `limit` | integer | no | Maximum number of results to return. Default `25` |

??? note "Upstream description"

    Search ENCODE chromatin accessibility experiments (ATAC-seq and DNase-seq). ATAC-seq and DNase-seq identify open chromatin regions including promoters, enhancers, and other regulatory elements. Returns experiment accessions, biosample details, and accessibility assay metadata. Open chromatin regions are sites of active gene regulation where transcription factors can bind. Use this to find chromatin accessibility data for specific tissues or cell types.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_search_chromatin_accessibility', 'arguments': {}},
)
```

## `ENCODE_search_hic_experiments` {#ENCODE_search_hic_experiments}

Searches ENCODE Hi-C experiments.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** hi-c, tad, chromatin loops, 3d genome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `assay_type` | string | no | Hi-C assay type: 'intact Hi-C' (newer protocol, better resolution) or 'Hi-C' (original protocol). Default is 'intact Hi-C'. Default `intact Hi-C` |
| `biosample_term_name` | string \| null | no | Biosample name filter (e.g., 'K562', 'GM12878', 'H1'). Leave empty for all. |
| `organism` | string | no | Organism scientific name (e.g., 'Homo sapiens', 'Mus musculus'). Default `Homo sapiens` |
| `limit` | integer | no | Maximum number of results to return. Default `25` |

??? note "Upstream description"

    Search ENCODE Hi-C and intact Hi-C experiments for 3D genome organization data. Hi-C captures chromatin interactions genome-wide, revealing topologically associating domains (TADs), chromatin loops, and A/B compartments. Returns experiment accessions, biosample details, and metadata. Use this to find 3D genome structure data for specific cell types or tissues. Intact Hi-C is the newer protocol with improved resolution.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_search_hic_experiments', 'arguments': {}},
)
```

## `ENCODE_search_histone_experiments` {#ENCODE_search_histone_experiments}

Searches ENCODE histone ChIP-seq experiments by mark and biosample.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** histone, h3k27ac, h3k4me3, chip-seq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `histone_mark` | string \| null | no | Histone modification mark to filter by (e.g., 'H3K4me3', 'H3K27ac', 'H3K27me3', 'H3K4me1', 'H3K36me3', 'H3K9me3'). Leave empty to search all histone marks. |
| `target` | string \| null | no | Alias for histone_mark. Histone modification mark (e.g., 'H3K27ac', 'H3K4me3'). |
| `biosample_term_name` | string \| null | no | Biosample name from ENCODE ontology (cell lines or tissues, NOT disease names). Examples: 'K562', 'HepG2', 'GM12878', 'liver', 'brain'. For disease-specific data, search by relevant cell lines (e.g., K562/MOLM-13 for AML). Leave empty to search all. |
| `biosample` | string \| null | no | Alias for biosample_term_name. Biosample (tissue or cell line, e.g., 'liver', 'K562'). |
| `organism` | string | no | Organism scientific name (e.g., 'Homo sapiens', 'Mus musculus'). Default `Homo sapiens` |
| `limit` | integer | no | Maximum number of results to return (1-100). Default `25` |
| `biosample_term` | string \| null | no | Alias for biosample_term_name. Biosample tissue or cell line (e.g., "breast epithelium", "liver", "K562"). |
| `cell_type` | string \| null | no | Alias for biosample_term_name. Cell type or tissue (e.g., 'GM12878', 'K562', 'liver', 'brain'). |
| `tissue` | string \| null | no | Alias for biosample_term_name. Tissue or cell line (e.g., 'liver', 'brain', 'breast epithelium', 'MCF-7'). |

??? note "Upstream description"

    Search ENCODE histone ChIP-seq experiments by histone modification mark, biosample, or organism. Returns experiment accessions, histone marks (H3K4me3, H3K27ac, H3K27me3, H3K36me3, H3K4me1, H3K9me3), biosample summaries, and metadata. Use this to find histone modification profiling data for specific marks or tissues. Common marks: H3K4me3 (active promoters), H3K27ac (active enhancers), H3K27me3 (polycomb repression), H3K4me1 (poised enhancers), H3K36me3 (gene bodies), H3K9me3 (heterochromatin).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_search_histone_experiments', 'arguments': {}},
)
```

## `ENCODE_search_methylation_experiments` {#ENCODE_search_methylation_experiments}

Searches ENCODE WGBS and RRBS methylation experiments.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** dna methylation, wgbs, rrbs, bisulfite
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `assay_type` | string | no | Methylation assay type: 'WGBS' (whole-genome bisulfite sequencing, comprehensive CpG coverage) or 'RRBS' (reduced-representation, enriched for CpG islands). Default `WGBS` |
| `biosample_term_name` | string \| null | no | Biosample name filter (e.g., 'K562', 'liver', 'brain', 'motor neuron'). Leave empty for all. |
| `organism` | string | no | Organism scientific name. Default `Homo sapiens` |
| `limit` | integer | no | Maximum number of results to return. Default `25` |

??? note "Upstream description"

    Search ENCODE whole-genome bisulfite sequencing (WGBS) and reduced-representation bisulfite sequencing (RRBS) experiments for DNA methylation profiling. Returns experiment accessions, biosample information, methylation assay type, and metadata. DNA methylation (5mC) at CpG sites is a key epigenetic mark regulating gene silencing, genomic imprinting, and X-inactivation. Use this to find methylation data for specific tissues, cell lines, or developmental stages.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_search_methylation_experiments', 'arguments': {}},
)
```

## `ENCORI_get_ceRNA_network` {#ENCORI_get_ceRNA_network}

Builds ENCORI ceRNA networks of genes and lncRNAs sharing miRNA binding with a query gene.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** cerna, mirna sponge, lncrna, encori
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string \| null | no | Query gene or lncRNA symbol whose ceRNA partners are wanted, e.g. 'PTEN', 'MALAT1' (aliases: 'ceRNA', 'gene_symbol'). |
| `assembly` | string \| null | no | Genome assembly (default 'hg38'). |
| `gene_type` | string \| null | no | Biotype of the query gene: 'mRNA' (default), 'lncRNA', 'pseudogene', 'sncRNA'. |
| `shared_mirna_min` | integer \| null | no | Minimum number of shared miRNA families required (default 5). |
| `pval` | number \| null | no | Maximum hypergeometric p-value (default 0.01). |
| `fdr` | number \| null | no | Maximum FDR (default 0.01). |
| `limit` | integer \| null | no | Maximum ceRNA partner rows to return (1-500, default 100). |

??? note "Upstream description"

    Build competing-endogenous-RNA (ceRNA / miRNA-sponge) networks via ENCORI (starBase) ceRNA module. Given a query gene/lncRNA (e.g. 'PTEN'), returns the genes and lncRNAs that share miRNA-binding families with it (candidate miRNA sponges), each with the number of shared miRNA families (shared_mirna_families), a p-value and an FDR. Lower pval/fdr and more shared families = stronger ceRNA crosstalk. Public, no authentication.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCORI_get_ceRNA_network', 'arguments': {}},
)
```

## `ENCORI_get_degradome_cleavage` {#ENCORI_get_degradome_cleavage}

Retrieves ENCORI degradome-seq validated miRNA cleavage sites by gene or miRNA.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** degradome, pare, mirna cleavage, encori
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string \| null | no | Target gene symbol to get cleaving miRNAs for, e.g. 'TP53' (alias: 'gene_symbol'). Mutually exclusive with 'mirna'. |
| `mirna` | string \| null | no | miRNA name to get cleavage targets for, e.g. 'hsa-miR-1-3p'. Mutually exclusive with 'gene'. |
| `assembly` | string \| null | no | Genome assembly (default 'hg19' — the only assembly with degradome data; hg38 returns nothing). |
| `gene_type` | string \| null | no | Target biotype: 'mRNA' (default), 'lncRNA', 'pseudogene', 'sncRNA'. |
| `degradome_exp_min` | integer \| null | no | Minimum number of supporting degradome experiments (default 1). |
| `clip_min` | integer \| null | no | Minimum number of supporting CLIP experiments (default 1). |
| `cell_type` | string \| null | no | Restrict to a cell line/tissue (default 'all'). |
| `limit` | integer \| null | no | Maximum cleavage rows to return (1-500, default 100). |

??? note "Upstream description"

    Retrieve miRNA cleavage sites validated by degradome-seq (PARE / Degradome) via ENCORI (starBase) degradomeRNA module: experimentally confirmed slicer cleavage of a target by a miRNA. Provide 'gene' (e.g. 'TP53') to get miRNAs that cleave it, or 'mirna' for a miRNA's cleavage targets. Each row reports the cleavage-event count, number of degradome experiments, degradome sites, total reads and category (I-IV signal confidence). NOTE: degradome data is built only for assembly 'hg19' (the default here); hg38 returns nothing. Public, no authentication.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCORI_get_degradome_cleavage', 'arguments': {}},
)
```

## `ENCORI_get_miRNA_targets` {#ENCORI_get_miRNA_targets}

Looks up ENCORI CLIP-supported and predicted miRNA-target interactions by miRNA or gene.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** mirna targets, encori, clip-seq, microrna
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `mirna` | string \| null | no | miRNA name to get targets for, e.g. 'hsa-miR-21-5p', 'hsa-let-7a-5p'. Mutually exclusive with 'gene'. |
| `gene` | string \| null | no | Gene symbol to get targeting miRNAs for, e.g. 'TP53', 'KRAS' (alias: 'gene_symbol'). Mutually exclusive with 'mirna'. |
| `assembly` | string \| null | no | Genome assembly to query (default 'hg38'). Use 'hg19' for the older human build. |
| `clip_min` | integer \| null | no | Minimum number of supporting CLIP-seq experiments (default 1 = require experimental support; set 0 to include prediction-only). |
| `program_min` | integer \| null | no | Minimum number of prediction programs that must call the interaction (default 1). |
| `limit` | integer \| null | no | Maximum interactions to return (1-500, default 50). |

??? note "Upstream description"

    Look up miRNA-target interactions from ENCORI (starBase): CLIP-seq-supported and computationally predicted targets of a miRNA, or the miRNAs that target a gene. Provide 'mirna' (e.g. 'hsa-miR-21-5p') to get its target genes, OR 'gene' (e.g. 'TP53') to get targeting miRNAs. Each result reports the number of supporting CLIP experiments (experimental evidence) and which prediction programs (TargetScan, miRanda, PITA, miRmap, microT, PicTar, RNA22) call it. Fills the gap where ToolUniverse had no miRNA target-lookup tool. Public, no authentication.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCORI_get_miRNA_targets', 'arguments': {}},
)
```

## `ENCORI_get_RBP_disease` {#ENCORI_get_RBP_disease}

Retrieves ENCORI RBP binding sites on genes that overlap COSMIC somatic mutations by tissue or disease.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** encori, starbase, rbp, somatic mutation, clip-seq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene` | string \| null | no | Bound target gene symbol, e.g. 'MYC' (alias: 'gene_symbol'). |
| `rbp` | string \| null | no | RNA-binding protein symbol, e.g. 'ACIN1' (alias: 'RBP'). |
| `tissue` | string \| null | no | Tissue filter, e.g. 'breast', 'lung', 'liver'. |
| `disease` | string \| null | no | Disease keyword, e.g. 'carcinoma', 'adenocarcinoma'. |
| `assembly` | string \| null | no | Genome assembly (default 'hg38'). |
| `limit` | integer \| null | no | Maximum association rows to return (1-500, default 100). |

??? note "Upstream description"

    Retrieve RBP-disease associations via ENCORI (starBase) RBPDisease module: RBP binding sites that overlap somatic mutations from cancer (COSMIC) in a given tissue/disease. Provide any of 'gene' (the bound target, e.g. 'MYC'), 'rbp', 'tissue' (e.g. 'breast'), or 'disease' (e.g. 'carcinoma'). Each row reports the disease name(s), COSMIC mutation IDs (disease_cosmic_id), mutation count (cosmic_num), affected sample count (sample_num) and supporting CLIP experiments. Public, no authentication.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCORI_get_RBP_disease', 'arguments': {}},
)
```

## `ENCORI_get_RBP_targets` {#ENCORI_get_RBP_targets}

Looks up ENCORI CLIP-seq RBP binding sites on transcripts, by RBP or target gene.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** encori, rbp, clip-seq, rna binding, transcript
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rbp` | string \| null | no | RNA-binding protein symbol to get targets for, e.g. 'PTBP1', 'ELAVL1' (alias: 'RBP'). Mutually exclusive with 'gene'. |
| `gene` | string \| null | no | Gene/transcript symbol to get binding RBPs for, e.g. 'TP53', 'MYC' (alias: 'gene_symbol'). Mutually exclusive with 'rbp'. |
| `assembly` | string \| null | no | Genome assembly (default 'hg38'; use 'hg19' for the older build). |
| `gene_type` | string \| null | no | Transcript biotype to query: 'mRNA' (default), 'lncRNA', 'pseudogene', 'sncRNA', 'miRNA'. |
| `clip_min` | integer \| null | no | Minimum number of supporting CLIP-seq experiments (default 1). |
| `pancancer_min` | integer \| null | no | Minimum pan-cancer support count (default 0). |
| `cell_type` | string \| null | no | Restrict to a CLIP cell line/tissue (default 'all'). |
| `limit` | integer \| null | no | Maximum binding rows to return (1-500, default 100). |

??? note "Upstream description"

    Look up RBP-RNA binding sites from CLIP-seq via ENCORI (starBase) RBPTarget module. Provide 'rbp' (e.g. 'PTBP1') to list the transcripts an RNA-binding protein binds, OR 'gene' (e.g. 'TP53') to list the RBPs that bind a transcript. Each result reports binding-cluster counts (cluster_num), the number of supporting CLIP experiments (total_clip_experiments) and CLIP sites (total_clip_sites), genomic location, and cell line/tissue. Public, no authentication.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCORI_get_RBP_targets', 'arguments': {}},
)
```

## `ENCORI_get_RNA_RNA_interactions` {#ENCORI_get_RNA_RNA_interactions}

Retrieves ENCORI RNA-RNA interaction partners of an RNA from PARIS/LIGR-seq/SPLASH data.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** encori, rna-rna interaction, lncrna, paris
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rna` | string \| null | no | RNA symbol whose duplex partners are wanted, e.g. 'MALAT1', 'NEAT1' (aliases: 'RNA', 'gene'). |
| `assembly` | string \| null | no | Genome assembly (default 'hg38'). |
| `gene_type` | string \| null | no | Biotype of the query RNA: 'lncRNA' (default), 'mRNA', 'pseudogene', 'sncRNA', 'miRNA'. Must match the RNA's real biotype or ENCORI returns no rows. |
| `interaction_min` | integer \| null | no | Minimum interaction count (default 1). |
| `exp_min` | integer \| null | no | Minimum number of supporting experiments (default 1). |
| `cell_type` | string \| null | no | Restrict to a cell line/tissue (default 'all'). |
| `limit` | integer \| null | no | Maximum interaction rows to return (1-500, default 100). |

??? note "Upstream description"

    Retrieve ncRNA-RNA (RNA-RNA duplex) interaction networks from high-throughput crosslinking (PARIS / LIGR-seq / SPLASH) via ENCORI (starBase) RNARNA module. Given an RNA (e.g. 'MALAT1'), returns its direct base-pairing partners with interaction/experiment/read support, predicted free energy and Smith-Waterman alignment score. NOTE: set 'gene_type' to the RNA's actual biotype (default 'lncRNA' suits MALAT1 and other lncRNAs); ENCORI silently returns no rows if the biotype is wrong. Public, no authentication.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCORI_get_RNA_RNA_interactions', 'arguments': {}},
)
```

## `ENCORI_scan_RBP_motifs` {#ENCORI_scan_RBP_motifs}

Scans ENCORI RBP binding motif enrichment by motif sequence or RBP.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** rbp motif, encori, clip peaks, motif scan
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `motif` | string \| null | no | RNA sequence motif to scan for (use U or T), e.g. 'UGCAUG'. Mutually exclusive with 'rbp'. |
| `rbp` | string \| null | no | RNA-binding protein symbol whose motifs are wanted, e.g. 'MBNL2' (alias: 'RBP'). |
| `assembly` | string \| null | no | Genome assembly (default 'hg38'). |
| `length` | string \| null | no | Motif length class: 'short' (default) or 'long'. |
| `rank_limit` | integer \| null | no | Per-dataset motif rank cutoff sent to ENCORI (1-100, default 10). |
| `limit` | integer \| null | no | Maximum motif rows to return (1-500, default 100). |

??? note "Upstream description"

    Scan RBP binding-motif enrichment via ENCORI (starBase) RBPMotifScan module. Provide 'motif' (e.g. 'UGCAUG') to find the RBPs and CLIP datasets whose peaks are enriched for that sequence, or 'rbp' to list an RBP's identified motifs. Each ranked row reports the identified motif, dataset ID, number of target peaks containing it (target_peak_num), percentage and enrichment p-value, plus a link to the HOMER motif matrix. Public, no authentication.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCORI_scan_RBP_motifs', 'arguments': {}},
)
```

## `ensembl_get_regulatory_features` {#ensembl_get_regulatory_features}

Gets Ensembl regulatory features such as promoters, enhancers and TF binding sites in a region.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `region`, `feature`
- **Keywords:** promoter, enhancer, tfbs, regulatory build
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `region` | string | yes | Genomic region in format 'chromosome:start-end' (e.g., '1:1000000-2000000') |
| `species` | string | no | Species name (default 'human') Default `human` |
| `feature` | string | yes | Feature type (must be 'regulatory' for this tool) |

??? note "Upstream description"

    Get regulatory features (promoters, enhancers, transcription factor binding sites) in a genomic region. Returns regulatory annotations including feature types, positions, and evidence. Use this to identify regulatory elements that may control gene expression.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ensembl_get_regulatory_features',
     'arguments': {'region': '...', 'feature': '...'}},
)
```

## `EnsemblReg_get_binding_matrix` {#EnsemblReg_get_binding_matrix}

Gets an Ensembl TF binding matrix by stable ID with frequencies and consensus.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `binding_matrix_id`
- **Keywords:** pwm, tf binding, motif, consensus
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | no | Species name. Use 'homo_sapiens' for human. Default: 'homo_sapiens'. Default `homo_sapiens` |
| `binding_matrix_id` | string | yes | Ensembl binding matrix stable ID. Format: 'ENSPFM0XXX'. Examples: 'ENSPFM0320' (PAX5), 'ENSPFM0423' (NR2F1), 'ENSPFM0001' (ALX3). |

??? note "Upstream description"

    Get a transcription factor binding matrix (position weight matrix/PWM) by stable ID from the Ensembl REST API. Returns the nucleotide frequency matrix, threshold, associated TF complexes, and consensus sequence. Binding matrices describe the DNA sequence preferences of transcription factors and are essential for predicting TF binding sites and understanding gene regulation. Example: ENSPFM0320 returns the PAX5 binding matrix with 18 positions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblReg_get_binding_matrix',
     'arguments': {'binding_matrix_id': '...'}},
)
```

## `EnsemblReg_get_motif_features` {#EnsemblReg_get_motif_features}

Gets TF binding motif features in a genomic region from Ensembl.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `region`
- **Keywords:** tf motif, binding site, genomic region
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | no | Species name. Use 'homo_sapiens' for human, 'mus_musculus' for mouse. Default: 'homo_sapiens'. Default `homo_sapiens` |
| `region` | string | yes | Genomic region in format 'chromosome:start-end'. Example: '7:140424943-140524564' for the BRAF locus, '17:7661779-7687538' for the TP53 locus. |

??? note "Upstream description"

    Get transcription factor (TF) binding motif features in a genomic region from the Ensembl REST API. Returns motif instances with their binding matrix IDs, TF complex names, scores, and genomic positions. Use this to identify which transcription factors may bind to a specific genomic region, essential for understanding gene regulation, epigenomics, and cis-regulatory elements. Example: BRAF locus (chr7:140424943-140524564) returns 50 TF binding motif instances including PAX5, NR2F1, CEBPA.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblReg_get_motif_features',
     'arguments': {'region': '...'}},
)
```

## `EnsemblReg_get_regulatory_elements` {#EnsemblReg_get_regulatory_elements}

Gets Ensembl Regulatory Build features (enhancers, promoters, CTCF sites, open chromatin) for a region.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `chrom`, `start`, `end`
- **Keywords:** regulatory build, enhancer, ctcf, open chromatin
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | no | Species name (e.g., 'homo_sapiens', 'mus_musculus'). Default `homo_sapiens` |
| `chrom` | string | yes | Chromosome number without 'chr' prefix (e.g., '17', '7', 'X'). |
| `start` | integer | yes | Start position (1-based). |
| `end` | integer | yes | End position (1-based). |

??? note "Upstream description"

    Get Ensembl regulatory features (enhancers, promoters, CTCF binding sites, open chromatin, TF binding sites) for a genomic region from the Ensembl Regulatory Build. The Ensembl Regulatory Build integrates epigenomic data from ENCODE, Roadmap Epigenomics, and BLUEPRINT to annotate regulatory elements across the human genome. Returns feature IDs, types, coordinates, and extended bounds. Use this to identify known regulatory elements in a region of interest for epigenomic analysis or variant interpretation.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EnsemblReg_get_regulatory_elements',
     'arguments': {'chrom': '...', 'start': 0, 'end': 0}},
)
```

## `EWASCatalog_search_by_cpg` {#EWASCatalog_search_by_cpg}

Gets EWAS Catalog associations for a CpG methylation site.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `cpg_id`
- **Keywords:** ewas, cpg, dna methylation, trait association
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `cpg_id` | string | yes | Illumina CpG probe identifier, e.g. 'cg05575921'. |
| `limit` | integer \| null | no | Max associations to return, 1-200. Default 50. |

??? note "Upstream description"

    Look up all published epigenome-wide association study (EWAS) results for one CpG methylation site from the EWAS Catalog (MRC-IEU). Returns trait, effect size, tissue, cohort, and publication for every association, ranked by significance. Example: cpg_id='cg05575921' (AHRR, the strongest known smoking-methylation signal).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EWASCatalog_search_by_cpg',
     'arguments': {'cpg_id': '...'}},
)
```

## `EWASCatalog_search_by_gene` {#EWASCatalog_search_by_gene}

Gets EWAS Catalog associations across CpG sites mapped to a gene.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `gene_symbol`
- **Keywords:** ewas, cpg, dna methylation, gene
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_symbol` | string | yes | HGNC gene symbol, e.g. 'AHRR', 'TP53'. |
| `limit` | integer \| null | no | Max associations to return, 1-200. Default 50. |

??? note "Upstream description"

    Look up all published epigenome-wide association study (EWAS) results annotated to one gene from the EWAS Catalog (MRC-IEU), across every CpG site mapped to it. Returns trait, effect size, tissue, cohort, and publication for each, ranked by significance. Example: gene_symbol='AHRR' (the smoking-methylation gene).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EWASCatalog_search_by_gene',
     'arguments': {'gene_symbol': '...'}},
)
```

## `FourDN_get_file_metadata` {#FourDN_get_file_metadata}

Gets 4DN file metadata for Hi-C contact matrices, TAD calls, A/B compartments and chromatin loops.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `file_accession`
- **Keywords:** 4dn, tad, compartments, chromatin loops, contact matrix
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | One of `get_file_metadata` Default `get_file_metadata` |
| `file_accession` | string | yes | 4DN file accession (e.g., '4DNFIIA7E3HL'). Obtain by searching with FourDN_search_data (item_type='File'). |
| `include_full_metadata` | boolean | no | Include complete API response in 'metadata' field (default: false). Set to true only if you need all raw fields. Omitting this reduces response size by ~97%. Default `False` |

??? note "Upstream description"

    Get detailed metadata for specific 4DN files including Hi-C contact matrices, TAD (Topologically Associating Domain) calls, A/B compartment annotations, and chromatin loop lists. Returns file type, format, size (check before downloading - files can be GBs), description, associated experiment, and biosource information. Prerequisites: Requires file accession - obtain by using FourDN_search_data with item_type='File'. Use for: checking file details before download, understanding processing pipeline, accessing quality metrics.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FourDN_get_file_metadata',
     'arguments': {'file_accession': '...'}},
)
```

## `FourDN_search_data` {#FourDN_search_data}

Searches the 4DN Data Portal for Hi-C chromatin conformation files and experiments by assay, cell type or keyword.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** 4dn, hi-c, micro-c, chromatin conformation, 3d genome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | One of `search` Default `search` |
| `query` | string | no | Search query (default: '*' for all) Default `*` |
| `item_type` | string | no | Type of item to search One of `File`, `Experiment`, `Biosource` Default `File` |
| `file_type` | string | no | Filter by file type. Common types: 'contact list' (processed contact matrices), 'pairs' (aligned read pairs), 'cooler' (multi-resolution contact matrices), 'mcool' (multi-resolution cooler), 'hic' (Juicer format). |
| `assay_title` | string | no | Filter by assay (e.g., 'Hi-C', 'Micro-C', 'Capture Hi-C') |
| `biosource_name` | string | no | Filter by cell type (e.g., 'H1-hESC', 'HFFc6', 'GM12878') |
| `limit` | integer | no | Maximum results to return Default `25` |

??? note "Upstream description"

    Search 4DN Data Portal for Hi-C chromatin conformation files and experiments. Access 350+ uniformly processed Hi-C contact files from multiple technologies (in situ Hi-C, Capture Hi-C, Micro-C XL). Search by assay type, cell type, or keywords. Returns file metadata with download URLs. Note: Downloads require free 4DN account (create at data.4dnucleome.org). Use for: finding Hi-C datasets, discovering 3D genome data, accessing TAD (Topologically Associating Domain) boundaries, A/B compartments, chromatin loop calls.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FourDN_search_data', 'arguments': {}},
)
```

## `GEO_search_atacseq_datasets` {#GEO_search_atacseq_datasets}

Searches NCBI GEO for ATAC-seq datasets.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** geo, atac-seq, chromatin accessibility
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string \| null | no | Search query (e.g., 'T cell differentiation', 'liver hepatocyte', 'hematopoiesis'). Leave empty to search broadly. |
| `organism` | string | no | Organism scientific name (e.g., 'Homo sapiens', 'Mus musculus'). Default `Homo sapiens` |
| `limit` | integer | no | Maximum number of results to return (1-100). Default `20` |

??? note "Upstream description"

    Search NCBI GEO for ATAC-seq chromatin accessibility datasets by query keyword, tissue, or cell type. Returns dataset accessions (GSE IDs), titles, summaries, organism, and sample counts. ATAC-seq measures open chromatin regions and is used to map regulatory elements, enhancers, and TF binding sites genome-wide.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GEO_search_atacseq_datasets', 'arguments': {}},
)
```

## `GEO_search_chipseq_datasets` {#GEO_search_chipseq_datasets}

Searches NCBI GEO for ChIP-seq datasets.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `query`
- **Keywords:** geo, chip-seq, tf binding, histone
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search terms for ChIP-seq datasets (e.g., 'H3K27ac liver', 'CTCF cancer', 'p53 ChIP-seq'). |
| `organism` | string | no | Organism filter. Default `Homo sapiens` |
| `limit` | integer | no | Maximum number of dataset IDs to return. Default `20` |

??? note "Upstream description"

    Search NCBI GEO for ChIP-seq (Chromatin Immunoprecipitation followed by sequencing) datasets. ChIP-seq profiles genome-wide binding of proteins (transcription factors, histones) to DNA. Returns dataset accessions, titles, summaries, and sample counts. Use this to find ChIP-seq studies for specific transcription factors, histone marks, or in specific cell types/diseases.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GEO_search_chipseq_datasets',
     'arguments': {'query': '...'}},
)
```

## `GEO_search_methylation_datasets` {#GEO_search_methylation_datasets}

Searches NCBI GEO for DNA methylation array datasets (450K, EPIC).

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `query`
- **Keywords:** geo, methylation array, 450k, epic
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search terms for methylation datasets (e.g., 'breast cancer methylation', 'brain methylation aging', 'CpG island methylation'). Combined with methylation platform filter. |
| `organism` | string | no | Organism filter (e.g., 'Homo sapiens', 'Mus musculus'). Default `Homo sapiens` |
| `limit` | integer | no | Maximum number of dataset IDs to return. Default `20` |

??? note "Upstream description"

    Search NCBI GEO for DNA methylation array datasets, including Illumina 450K, EPIC (850K), and other methylation profiling platforms. Returns dataset accessions (GSE IDs), titles, summaries, platform information, sample counts, and organisms. DNA methylation arrays measure methylation levels at hundreds of thousands of CpG sites genome-wide. Use this to find published methylation studies for specific conditions, tissues, or diseases. Follow up with GEO_get_dataset_details for full metadata.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GEO_search_methylation_datasets',
     'arguments': {'query': '...'}},
)
```

## `HOCOMOCO_get_motif` {#HOCOMOCO_get_motif}

Gets a HOCOMOCO TF binding motif with gene names, TF class, quality, consensus and optional PWM.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `motif_id`
- **Keywords:** tf motif, pwm, transcription factor, hocomoco, consensus
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `motif_id` | string | yes | HOCOMOCO motif identifier. Format: GENE.COLLECTION.SUBTYPE.SOURCE.QUALITY. Examples: 'P53.H14CORE.0.P.B', 'CTCF.H14CORE.0.P.B', 'STAT3.H14CORE.0.P.B'. |
| `include_pwm` | boolean | no | Whether to include the position weight matrix (PWM) in the response. Default: false. Default `False` |

??? note "Upstream description"

    Get detailed information about a specific HOCOMOCO transcription factor binding motif. Returns gene names, UniProt IDs, TF classification, quality metrics, consensus sequence, and optionally the position weight matrix (PWM). Quality grades: A (best), B, C, D (lowest). Data sources: P (ChIP-Seq peaks), S (HT-SELEX). Example: 'P53.H14CORE.0.P.B' returns TP53 motif details.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HOCOMOCO_get_motif', 'arguments': {'motif_id': '...'}},
)
```

## `HOCOMOCO_search_motifs` {#HOCOMOCO_search_motifs}

Searches HOCOMOCO v14 TF binding motifs by gene/protein name returning IDs, quality and consensus.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `query`
- **Keywords:** tf motif, transcription factor binding, chip-seq, hocomoco
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Gene symbol or protein name to search for TF binding motifs. Examples: 'TP53', 'CTCF', 'STAT3', 'SOX2', 'MYC'. |

??? note "Upstream description"

    Search HOCOMOCO v14 for transcription factor binding motifs by gene or protein name. Returns motif IDs, quality scores, consensus sequences, and TF classification. HOCOMOCO contains high-quality TF binding models derived from ChIP-Seq data for human and mouse. Complementary to JASPAR. Example: query 'TP53' returns P53.H14CORE.0.P.B with consensus sequence and quality grade.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'HOCOMOCO_search_motifs', 'arguments': {'query': '...'}},
)
```

## `jaspar_get_matrix` {#jaspar_get_matrix}

Returns metadata and the PFM/PWM for a single JASPAR matrix by matrix ID.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `matrix_id`
- **Keywords:** pfm, pwm, tf motif, binding site
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `matrix_id` | string | yes | JASPAR matrix ID (e.g., 'MA0002.2'). You can discover IDs with `JASPAR_get_transcription_factors`. |

??? note "Upstream description"

    Get full details for a single JASPAR matrix by `matrix_id` (e.g., MA0002.2). Returns metadata plus PFM/PWM. For reproducible analyses, use `jaspar_get_matrix_versions` with the corresponding `base_id` (e.g., MA0002) to see all versions and choose/cite the exact version you used.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'jaspar_get_matrix', 'arguments': {'matrix_id': '...'}},
)
```

## `jaspar_get_matrix_versions` {#jaspar_get_matrix_versions}

Lists all versions of a JASPAR TF binding matrix for a base ID.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `base_id`
- **Keywords:** motif version, matrix base id, tf binding, reproducibility
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `base_id` | string | yes | JASPAR base matrix ID without version (e.g., 'MA0002'). You can get `base_id` from `jaspar_get_matrix`. |
| `page` | integer | no | Page number (1-based). Default `1` |
| `page_size` | integer | no | Results per page. Default `20` |

??? note "Upstream description"

    List all versions available for a JASPAR matrix base ID (e.g., base_id MA0002). This is useful for reproducibility (choose a specific matrix version like MA0002.2) and for tracking updates across releases.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'jaspar_get_matrix_versions',
     'arguments': {'base_id': '...'}},
)
```

## `JASPAR_get_transcription_factors` {#JASPAR_get_transcription_factors}

Lists all JASPAR transcription factor binding profile matrices with pagination and optional collection filter.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** transcription factor, binding profile, motif, pfm
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `collection` | string | no | JASPAR collection (e.g., CORE). Default `CORE` |
| `page` | integer | no | Page number (1-based). Default `1` |
| `page_size` | integer | no | Results per page (JASPAR `page_size`). Default `20` |

??? note "Upstream description"

    List all JASPAR transcription factor binding profiles. Returns all matrices; does NOT support filtering by TF name. To search by name, use jaspar_search_matrices with \{"name": "FOXA1"\}. Supports pagination via page/page_size.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'JASPAR_get_transcription_factors', 'arguments': {}},
)
```

## `jaspar_search_matrices` {#jaspar_search_matrices}

Searches and filters JASPAR TF binding matrices by text, TF name, collection, taxon or species.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** tf motif search, transcription factor, pwm, binding site
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `search` | string | no | Free-text search across matrices (e.g., 'CTCF'). |
| `name` | string | no | Filter by TF name (e.g., 'CTCF'). |
| `collection` | string | no | Filter by collection (e.g., 'CORE'). |
| `tax_group` | string | no | Filter by taxonomic group (e.g., 'vertebrates'). |
| `species` | string | no | Filter by NCBI taxonomy ID (e.g., '9606' for human). Sent to JASPAR as its `tax_id` query parameter. |
| `page` | integer | no | Page number (1-based). Default `1` |
| `page_size` | integer | no | Results per page. Default `20` |

??? note "Upstream description"

    Search/filter JASPAR matrices using common query parameters. Typical filters include `search` (free text), `name`, `collection`, `tax_group`, `species`, and paging via `page`/`page_size`.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'jaspar_search_matrices', 'arguments': {}},
)
```

## `MEME_discover_motifs` {#MEME_discover_motifs}

Runs MEME de novo motif discovery on a set of DNA sequences.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `sequences`
- **Keywords:** de novo motif, meme, tf binding site, fasta
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `discover_motifs` Default `discover_motifs` |
| `sequences` | string | yes | Multiple DNA sequences in FASTA format (minimum 2 sequences). Should contain related sequences expected to share common motifs, e.g., promoter regions of co-regulated genes. |
| `nmotifs` | integer \| null | no | Maximum number of motifs to find. Default 3. Default `3` |
| `minw` | integer \| null | no | Minimum motif width in nucleotides. Default 6. Default `6` |
| `maxw` | integer \| null | no | Maximum motif width in nucleotides. Default 50. Default `50` |
| `distribution` | string \| null | no | Motif site distribution model. 'oops' = one per sequence, 'zoops' = zero or one per sequence (default), 'anr' = any number of repetitions. One of `oops`, `zoops`, `anr` Default `zoops` |
| `scan_rc` | boolean \| null | no | If true (default), search both strands. Default `True` |

??? note "Upstream description"

    Run de novo motif discovery on a set of DNA sequences using the MEME algorithm from the MEME Suite. MEME finds overrepresented ungapped sequence patterns (motifs) that may correspond to transcription factor binding sites, splice sites, or other functional elements. Requires at least 2 FASTA sequences as input. Returns discovered motifs with consensus sequences, E-values, probability matrices, and site counts. Note: this submits a job to meme-suite.org and polls for results (typically 30-120 seconds depending on input size).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MEME_discover_motifs', 'arguments': {'sequences': '...'}},
)
```

## `MEME_fimo_scan` {#MEME_fimo_scan}

Scans DNA sequences for occurrences of a MEME-format motif using FIMO and returns significant sites.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `sequences`, `motif_text`
- **Keywords:** fimo, motif scan, tf binding, p-value
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `fimo_scan` Default `fimo_scan` |
| `sequences` | string | yes | DNA sequences in FASTA format to scan for motif occurrences. Example: '&gt;seq1\\nATGCATGCTATAAAAGGATCC\\n&gt;seq2\\nGGATCCATGCATGCATGCAT' |
| `motif_text` | string | yes | Motif(s) in MEME format. Must include the full MEME header (MEME version 5, ALPHABET, strands, Background letter frequencies) followed by MOTIF blocks with letter-probability matrices. Example for a TATA box motif: 'MEME version 5\\n\\nALPHABET= ACGT\\n\\nstrands: + -\\n\\nBackground letter frequencies\\nA 0.25 C 0.25 G 0.25 T 0.25\\n\\nMOTIF TATABOX\\nletter-probability matrix: alength= 4 w= 6 nsites= 100… |
| `pvalue_threshold` | number \| null | no | P-value threshold for reporting motif occurrences. Default 1e-4. Lower values are more stringent. Default `0.0001` |
| `scan_rc` | boolean \| null | no | If true (default), scan both strands. Set false to scan only the given strand. Default `True` |

??? note "Upstream description"

    Scan DNA sequences for transcription factor (TF) binding motifs using FIMO (Find Individual Motif Occurrences) from the MEME Suite. Provide sequences in FASTA format and a motif in MEME format, and FIMO returns all statistically significant binding sites with p-values, scores, strand, and matched sequences. Ideal for scanning promoter regions for known TF binding sites from JASPAR or HOCOMOCO. Note: this submits a job to meme-suite.org and polls for results (typically 10-30 seconds).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MEME_fimo_scan',
     'arguments': {'sequences': '...', 'motif_text': '...'}},
)
```

## `MEME_tomtom_compare` {#MEME_tomtom_compare}

Compares a query motif against known TF motif databases with TOMTOM and returns similar motifs.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `query_motif`
- **Keywords:** tomtom, motif comparison, tf annotation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `tomtom_compare` Default `tomtom_compare` |
| `query_motif` | string | yes | Query motif in MEME format. Must include the full MEME header followed by one or more MOTIF blocks with letter-probability matrices. |
| `target_db` | string \| null | no | Target motif database to search against. Options: 'JASPAR2026_vertebrates' (default, recommended), 'JASPAR2026_all', 'JASPAR2024_vertebrates', 'HOCOMOCO_v12'. Use MEME_list_databases to see all options. Default `JASPAR2026_vertebrates` |
| `evalue_threshold` | number \| null | no | E-value threshold for reporting matches. Default 0.5. Must be &lt;= 1.0. Default `0.5` |
| `comparison_function` | string \| null | no | Distance function for motif comparison. 'pearson' = Pearson correlation (default), 'ed' = Euclidean distance, 'sandelin' = Sandelin-Wasserman, 'allr' = average log-likelihood ratio. One of `pearson`, `ed`, `sandelin`, `allr`, `allr5` Default `pearson` |

??? note "Upstream description"

    Compare a query motif against known motif databases (JASPAR, HOCOMOCO, CIS-BP) using TOMTOM from the MEME Suite. Given a motif in MEME format, identifies the most similar known transcription factor binding motifs. Useful for annotating de novo discovered motifs — e.g., after running MEME motif discovery, use TOMTOM to determine which TF likely binds the discovered motif. Returns matches ranked by E-value with query/target consensus sequences. Note: this submits a job to meme-suite.org and polls for results (typically 5-20 seconds).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MEME_tomtom_compare',
     'arguments': {'query_motif': '...'}},
)
```

## `OmniPath_get_dorothea_regulon` {#OmniPath_get_dorothea_regulon}

Returns the DoRothEA regulon of a transcription factor with target genes, mode of regulation and confidence.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `tf_gene`
- **Keywords:** dorothea, regulon, tf targets, mode of regulation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `tf_gene` | string | yes | Transcription factor gene symbol (e.g., 'TP53', 'MYC', 'STAT3', 'NF2') |
| `confidence_levels` | string \| null | no | Comma-separated DoRothEA confidence levels to include (e.g., 'A,B' for high confidence only). Options: A, B, C, D, E. Default: all levels. |

??? note "Upstream description"

    Get the complete DoRothEA regulon for a transcription factor from OmniPath. Returns all target genes with mode of regulation (MoR: +1=activation, -1=repression), DoRothEA confidence levels, and source databases. Optionally filter by confidence levels (A, B, C). DoRothEA is the standard resource for TF activity estimation in tools like VIPER, decoupleR, and PROGENy. Use for: TF regulon retrieval for activity scoring, network analysis.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OmniPath_get_dorothea_regulon',
     'arguments': {'tf_gene': '...'}},
)
```

## `OmniPath_get_tf_target_interactions` {#OmniPath_get_tf_target_interactions}

Returns TF-target gene interactions from DoRothEA and CollecTRI with mode of regulation and confidence.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** tf-target, collectri, dorothea, gene regulation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `tf_gene` | string \| null | no | Transcription factor gene symbol (e.g., 'TP53', 'MYC', 'STAT3'). Mutually exclusive with target_gene or can be combined. |
| `target_gene` | string \| null | no | Target gene symbol to find regulating TFs (e.g., 'CDKN1A', 'BCL2'). Mutually exclusive with tf_gene or can be combined. |
| `confidence_level` | string \| null | no | DoRothEA confidence level filter: A=highest (1-2 sources), B=high (3+ sources), C=medium, D=low, E=lowest. A/B recommended for analysis. One of `A`, `B`, `C`, `D`, `E`, `None` |

??? note "Upstream description"

    Get transcription factor (TF) to target gene interactions from OmniPath DoRothEA and CollecTRI regulons. Returns directed TF-target pairs with mode of regulation (activation/repression), DoRothEA confidence levels (A-E), and supporting evidence. Query by TF gene symbol or target gene. Use for: TF activity inference, regulatory network construction, identifying TF targets for VIPER/decoupleR analysis.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OmniPath_get_tf_target_interactions', 'arguments': {}},
)
```

## `RegulomeDB_query_variant` {#RegulomeDB_query_variant}

Returns RegulomeDB regulatory evidence and score for a variant by rsID.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** regulomedb score, regulatory variant, rsid, tf binding
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `rsid` | string | no | dbSNP rsID (e.g., 'rs4994', 'rs429358') |
| `variant` | string | no | Alias for rsid. dbSNP rsID (e.g., 'rs7903146', 'rs429358'). |

??? note "Upstream description"

    Query RegulomeDB for regulatory variant annotations using rsID

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'RegulomeDB_query_variant', 'arguments': {}},
)
```

## `ReMap_get_peaks_in_region` {#ReMap_get_peaks_in_region}

Returns ReMap ChIP-seq transcriptional regulator peaks overlapping a genomic interval.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `region`
- **Keywords:** chip-seq, tf binding, peaks, genomic interval
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation selector. Must be 'get_peaks_in_region'. One of `get_peaks_in_region` Default `get_peaks_in_region` |
| `region` | string | yes | Genomic interval as chrom:start-end (e.g. 'chr1:1000000-1100000'). Commas in coordinates are tolerated. |
| `assembly` | string | no | Genome assembly. One of `hg38`, `hg19`, `mm10`, `mm39` Default `hg38` |
| `version` | string | no | ReMap catalog release year. One of `2022`, `2020`, `2018` Default `2022` |
| `datatype` | string | no | Peak datatype: 'all' for the non-redundant merged peak set across all TFs. Default `all` |
| `limit` | integer \| null | no | Optional cap on the number of peaks returned (peak_count still reports the true total). Omit to return all overlapping peaks. |

??? note "Upstream description"

    Retrieve all ChIP-seq transcriptional-regulator (TR) binding peaks from the ReMap catalog that overlap a genomic interval. Give a region (chrom:start-end) and get back every overlapping peak with its genomic coordinates plus the bound TF, the cell line/tissue (biotype), the source experiment accession (GEO/ENCODE), and any treatments. This is the defining ReMap capability (region in -&gt; all overlapping TR peaks out). Use for: enumerating which transcription factors bind a locus, finding regulatory occupancy of a promoter/enhancer interval, listing TFs near a variant. Keep regions modest (&lt;=1 Mb) - large intervals can return tens of thousands of peaks.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ReMap_get_peaks_in_region',
     'arguments': {'region': '...'}},
)
```

## `ReMap_get_transcription_factor_binding` {#ReMap_get_transcription_factor_binding}

Returns ReMap transcription factor binding sites for a gene and cell type.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `gene_name`
- **Keywords:** tfbs, chip-seq, cell type
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_name` | string | yes | Gene symbol (e.g., BRCA1, TP53, MYC) |
| `cell_type` | string | no | Cell type (e.g., HepG2, K562, MCF-7) Default `HepG2` |
| `limit` | integer | no | Number of results to return Default `10` |

??? note "Upstream description"

    Get transcription factor binding sites from ReMap database for specific genes and cell types

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ReMap_get_transcription_factor_binding',
     'arguments': {'gene_name': '...'}},
)
```

## `ReMap_list_datasets_for_target` {#ReMap_list_datasets_for_target}

Lists ReMap ChIP-seq datasets for a transcription factor with optional peak counts.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `target`
- **Keywords:** chip-seq, geo, encode, dataset
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `target` | string | yes | Transcription factor / target name, e.g. 'FOXA1'. |
| `taxid` | integer | no | NCBI taxonomy id. Default 9606 (human). Default `9606` |
| `experiment` | string | no | Optional experiment accession prefix to filter on, e.g. 'GSE23852'. Dataset names look like 'GSE23852.FOXA1.MCF-7_ETOH'. |
| `count_peaks` | boolean | no | Download each dataset's BED and report its peak count. Slower; use when the question asks how many peaks. Default `False` |

??? note "Upstream description"

    List ReMap ChIP-seq datasets for a transcription factor, one row per dataset, optionally filtered to a GEO/ENCODE experiment accession and optionally with each dataset's peak count. A GEO series is NOT one dataset: GSE23852 for FOXA1 holds two (MCF-7_ETOH and MCF-7_E2) with different peak counts, so 'how many peaks in experiment X' usually needs the per-dataset split rather than a single number. Set count_peaks=true to download each dataset's BED and count peaks (slower).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ReMap_list_datasets_for_target',
     'arguments': {'target': '...'}},
)
```

## `SCREEN_get_cCRE_by_accession` {#SCREEN_get_cCRE_by_accession}

Looks up ENCODE SCREEN candidate cis-regulatory elements by accession and returns their classification (promoter-like, enhancer-like, CTCF) and genomic location.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `accession`
- **Keywords:** encode, enhancer-like, promoter-like, ctcf, regulatory element, grch38
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string \| array | yes | A single SCREEN cCRE accession string, or a list of accession strings (e.g. 'EH38E2666166'). |
| `assembly` | string | no | Genome assembly: 'GRCh38' (human, default) or 'mm10' (mouse). Default `GRCh38` |

??? note "Upstream description"

    Look up one or more candidate cis-Regulatory Elements (cCREs) in the ENCODE SCREEN registry by SCREEN accession (e.g. 'EH38E2666166'). Returns the element classification (PLS=promoter-like, pELS/dELS=proximal/distal enhancer-like, CTCF-only, DNase-H3K4me3), its genomic coordinates, the underlying representative DNase hypersensitive site (rDHS) accession, and whether it is CTCF-bound. Default assembly is GRCh38. Accessions are obtained from SCREEN_search_cCREs_by_region. Reliability note: queries the SCREEN staging GraphQL endpoint (the production host is rate-limited/unavailable for the API), so it may be intermittently unstable.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SCREEN_get_cCRE_by_accession',
     'arguments': {'accession': '...'}},
)
```

## `SCREEN_get_regulatory_elements` {#SCREEN_get_regulatory_elements}

Retrieves ENCODE SCREEN candidate cis-regulatory elements (enhancer, promoter or insulator types) associated with a gene symbol.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** ccre, enhancer, promoter, encode, cis-regulatory element, insulator
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_name` | string | no | Gene symbol to search for regulatory elements (e.g., BRCA1, TP53). Alias: 'query'. |
| `query` | string | no | Gene symbol to search for regulatory elements (alias for gene_name). |
| `element_type` | string | no | Type of regulatory element (promoter, enhancer, insulator) Default `enhancer` |
| `limit` | integer | no | Number of results to return Default `10` |

??? note "Upstream description"

    Get candidate cis-regulatory elements (cCREs) from SCREEN database for specific genomic regions

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SCREEN_get_regulatory_elements', 'arguments': {}},
)
```

## `SCREEN_search_cCREs_by_region` {#SCREEN_search_cCREs_by_region}

Queries the ENCODE SCREEN registry for candidate cis-regulatory elements overlapping a chromosome region, optionally filtered by element type.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `chrom`, `start`, `end`
- **Keywords:** encode, genomic region, enhancer, promoter, ccre registry, chromosome coordinates
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chrom` | string | yes | Chromosome, e.g. 'chr8' or '8' (the 'chr' prefix is added automatically). |
| `start` | integer | yes | Region start coordinate (1-based). |
| `end` | integer | yes | Region end coordinate (1-based, must be greater than start). |
| `assembly` | string | no | Genome assembly: 'GRCh38' (human, default) or 'mm10' (mouse). Default `GRCh38` |
| `element_type` | string \| null | no | Optional element-type hint to bias results toward a class: 'PLS', 'pELS', 'dELS', 'CTCF-only', or 'DNase-H3K4me3'. Omit (null) to return all overlapping cCREs. Default `None` |

??? note "Upstream description"

    Query the ENCODE SCREEN Registry of candidate cis-Regulatory Elements (cCREs) overlapping a genomic region. Returns experimentally grounded regulatory annotations: each cCRE's element type (PLS=promoter-like, pELS/dELS=proximal/distal enhancer-like, CTCF-only, DNase-H3K4me3) and its DNase, H3K4me3 (promoter), H3K27ac (enhancer), and CTCF z-scores aggregated across ENCODE biosamples. This is the curated regulatory-genomics counterpart to sequence-to-track deep-learning predictors (Enformer/Borzoi), which have no public REST endpoint. Use it to find which regulatory elements lie in a locus and how strong their chromatin signals are. Coordinates are genome-build specific; default assembly is GRCh38. Reliability note: queries the SCREEN staging GraphQL endpoint (the production host is rate-limited/unavailable for the API), so it may be intermittently unstable.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SCREEN_search_cCREs_by_region',
     'arguments': {'chrom': '...', 'start': 0, 'end': 0}},
)
```

## `SGD_get_regulation` {#SGD_get_regulation}

Gets the curated transcriptional regulation network for a yeast gene from SGD (regulator-target pairs with evidence and direction).

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `sgd_id`
- **Keywords:** yeast, transcription factor, regulatory network, regulator target, gene regulation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sgd_id` | string | yes | SGD locus identifier. Examples: 'S000000364' (CDC28), 'S000000259' (GAL4). |

??? note "Upstream description"

    Get the transcriptional regulation network for a yeast gene from SGD: curated regulator-&gt;target relationships with evidence, direction, and references. Returns pairs where the queried gene is regulator or target. Example: S000000364 (CDC28) is regulated by IXR1 (S000001515) with 'DNA to cDNA expression microarray evidence', regulation_of=transcription.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SGD_get_regulation', 'arguments': {'sgd_id': '...'}},
)
```

## `UCSC_get_cpg_islands` {#UCSC_get_cpg_islands}

Gets UCSC CpG island annotations for a genomic region.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `chrom`, `start`, `end`
- **Keywords:** cpg island, ucsc, promoter, methylation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string | no | Genome assembly (e.g., 'hg38', 'hg19', 'mm10', 'mm39'). Default `hg38` |
| `chrom` | string | yes | Chromosome name (e.g., 'chr17', 'chr1', 'chrX'). |
| `start` | integer | yes | Start position (0-based, inclusive). |
| `end` | integer | yes | End position (0-based, exclusive). |

??? note "Upstream description"

    Get CpG island annotations for a genomic region from UCSC Genome Browser. CpG islands are genomic regions with high CpG dinucleotide density, often found at gene promoters. Methylation of CpG islands is associated with gene silencing. Returns CpG island locations, length, CpG count, GC content, and observed/expected CpG ratio. Use this to identify CpG islands near genes of interest for methylation analysis, or to characterize the regulatory landscape of a genomic region.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UCSC_get_cpg_islands',
     'arguments': {'chrom': '...', 'start': 0, 'end': 0}},
)
```

## `UCSC_get_encode_cCREs` {#UCSC_get_encode_cCREs}

Gets ENCODE4 candidate cis-regulatory elements for a region from UCSC.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `chrom`, `start`, `end`
- **Keywords:** ccre, ucsc, enhancer-like, promoter-like
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string | no | Genome assembly (e.g., 'hg38', 'mm10'). Default `hg38` |
| `chrom` | string | yes | Chromosome name (e.g., 'chr17'). |
| `start` | integer | yes | Start position (0-based). |
| `end` | integer | yes | End position (0-based). |

??? note "Upstream description"

    Get ENCODE4 candidate cis-Regulatory Elements (cCREs) for a genomic region from UCSC Genome Browser. cCREs are classified by epigenetic signatures into: Promoter-like (PLS, red), Proximal enhancer-like (pELS, orange), Distal enhancer-like (dELS, yellow), DNase-H3K4me3 (blue), CTCF-only/CTCF-bound (green). Each cCRE includes max Z-scores for DNase, H3K4me3, H3K27ac, and CTCF signals across biosamples. Use this to identify regulatory elements in a region of interest for functional genomics or GWAS variant interpretation.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UCSC_get_encode_cCREs',
     'arguments': {'chrom': '...', 'start': 0, 'end': 0}},
)
```

## `UCSC_get_tf_binding_clusters` {#UCSC_get_tf_binding_clusters}

Gets ENCODE3 TF ChIP-seq clusters for a region from UCSC.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `chrom`, `start`, `end`
- **Keywords:** tf binding, chip-seq clusters, ucsc, encode3
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `genome` | string | no | Genome assembly (e.g., 'hg38', 'hg19'). Default `hg38` |
| `chrom` | string | yes | Chromosome name (e.g., 'chr17'). |
| `start` | integer | yes | Start position (0-based). |
| `end` | integer | yes | End position (0-based). |

??? note "Upstream description"

    Get Transcription Factor ChIP-seq Clusters from ENCODE3 for a genomic region via UCSC Genome Browser. Returns clusters of TF binding sites from 340 factors across 129 cell types, aggregated from ENCODE3 ChIP-seq peak data. Each cluster shows the transcription factor name, peak score, and number of source experiments. Use this to identify which transcription factors bind to a specific genomic region, discover regulatory hotspots, or interpret GWAS variants by finding overlapping TF binding sites.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UCSC_get_tf_binding_clusters',
     'arguments': {'chrom': '...', 'start': 0, 'end': 0}},
)
```

## `UniBind_get_dataset` {#UniBind_get_dataset}

Retrieves UniBind direct transcription factor-DNA binding-site details for one dataset, including TF, cell type and motif information.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** `dataset_id`
- **Keywords:** tf binding sites, chip-seq, transcription factor, motif, tfbs
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation selector. Must be 'get_dataset'. One of `get_dataset` Default `get_dataset` |
| `dataset_id` | string | yes | UniBind dataset identifier, e.g. 'EXP030726.neural_stem_cells.SMAD3'. Obtain from UniBind_search_datasets (the 'dataset_id' field). |

??? note "Upstream description"

    Retrieve full direct transcription factor (TF)-DNA binding-site detail for one UniBind dataset by its dataset_id (e.g. 'EXP030726.neural_stem_cells.SMAD3', obtained from UniBind_search_datasets). Returns the TF name, cell line/tissue, biological condition, source identifiers (GEO/ENCODE/PubMed), the JASPAR motif IDs used, and a list of TFBS model blocks. Each TFBS model gives the JASPAR motif (jaspar_id + version), total number of predicted binding sites (total_tfbs), the PWM score threshold, the centrality distance threshold, the CentriMo adjusted enrichment p-value, and direct download URLs for the binding sites in BED format, the underlying sequences in FASTA, and a summary plot. These are motif-anchored, base-pair-resolution, experimentally supported direct binding sites - distinct from raw ChIP-seq peaks. Public, no API key.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UniBind_get_dataset', 'arguments': {'dataset_id': '...'}},
)
```

## `UniBind_list_tfs` {#UniBind_list_tfs}

Lists transcription factors profiled in UniBind, optionally filtered by TF gene symbol substring.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** transcription factors, tf catalog, binding sites, chip-seq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation selector. Must be 'list_tfs'. One of `list_tfs` Default `list_tfs` |
| `search` | string \| null | no | Case-insensitive substring filter applied client-side to TF gene symbols, e.g. 'SMAD', 'FOX', 'GATA'. Optional; omit to list all TFs. |
| `limit` | integer | no | Maximum number of TF names to return after filtering. Default `200` |

??? note "Upstream description"

    List the catalog of transcription factors (TFs) profiled in the UniBind database of direct TF-DNA binding sites. Optionally filter by a case-insensitive substring of the TF gene symbol via 'search' (e.g. 'SMAD' returns SMAD2, SMAD3, SMAD4...). Returns a sorted, de-duplicated list of TF gene symbols. Use this to discover valid tf_name values for UniBind_search_datasets. Public, no API key.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UniBind_list_tfs', 'arguments': {}},
)
```

## `UniBind_search_datasets` {#UniBind_search_datasets}

Searches UniBind TF-DNA binding-site datasets (one ChIP-seq experiment per dataset) by TF symbol, cell type or species.

- **Category:** `regulation_epigenome` (Regulatory elements, TF binding, chromatin, methylation and 3D genome)
- **Required arguments:** none
- **Keywords:** tf binding, chip-seq, transcription factor, cell type, tfbs
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation selector. Must be 'search_datasets'. One of `search_datasets` Default `search_datasets` |
| `tf_name` | string \| null | no | Transcription factor gene symbol to filter by, e.g. 'CTCF', 'SMAD3', 'ESR1'. Exact match on UniBind's TF name. Optional. |
| `species` | string \| null | no | Species scientific name to filter by, e.g. 'Homo sapiens', 'Mus musculus', 'Drosophila melanogaster'. Optional. |
| `cell_line` | string \| null | no | Cell type / tissue / condition name to filter by, e.g. 'neural stem cells', 'MCF7'. Optional. |
| `collection` | string \| null | no | UniBind collection: 'Robust' (high-confidence, metadata-curated) or 'Permissive' (broader). Optional. One of `Robust`, `Permissive`, `None` |
| `order` | string \| null | no | Field to order results by, e.g. 'tf_name' or '-tf_name' (prefix '-' for descending). Optional. |
| `page` | integer | no | 1-based page number for pagination. Default `1` |
| `page_size` | integer | no | Number of datasets per page (max 1000). Default `25` |

??? note "Upstream description"

    Search the UniBind catalog of curated, experimentally derived direct transcription factor (TF)-DNA binding-site datasets (one dataset = one ChIP-seq experiment for one TF in one cell type/condition). Filter by TF gene symbol (tf_name), species scientific name (species, e.g. 'Homo sapiens'), cell type/tissue (cell_line), and collection ('Robust' or 'Permissive'). All filters are optional and compose. Returns a paginated list of matching datasets with their TF name, total peak count, and dataset_id (use UniBind_get_dataset for full binding-site detail). UniBind provides motif-anchored, base-pair-resolution direct binding sites derived from ChIP-seq via the DAMO/ChIP-eat pipeline with JASPAR motifs - distinct from JASPAR/HOCOMOCO motif matrices and from ReMap raw ChIP-seq peaks. Public, no API key.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UniBind_search_datasets', 'arguments': {}},
)
```
