# ToolUniverse: sequencing_data (`tooluniverse`)

DNAKit's MCP tool index keeps **845** DNA/gene tools of this server (2927 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

This page lists the 91 tools in category `sequencing_data`; other categories: [`sequence`](tooluniverse-sequence.md) (44), [`genome_annotation`](tooluniverse-genome_annotation.md) (253), [`variant`](tooluniverse-variant.md) (155), [`association`](tooluniverse-association.md) (138), [`expression`](tooluniverse-expression.md) (56), [`regulation_epigenome`](tooluniverse-regulation_epigenome.md) (54), [`prediction_model`](tooluniverse-prediction_model.md) (17), [`molecular_tools`](tooluniverse-molecular_tools.md) (23), [`breeding_population`](tooluniverse-breeding_population.md) (14).

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
| [`arrayexpress_get_experiment`](#arrayexpress_get_experiment) | `sequencing_data` | Retrieves ArrayExpress experiment metadata by accession (E-MTAB, E-GEOD). |
| [`arrayexpress_get_experiment_files`](#arrayexpress_get_experiment_files) | `sequencing_data` | Lists raw, processed and MAGE-TAB files with download paths for an ArrayExpress experiment. |
| [`arrayexpress_get_experiment_samples`](#arrayexpress_get_experiment_samples) | `sequencing_data` | Returns sample characteristics and experimental factors for an ArrayExpress experiment. |
| [`arrayexpress_search_experiments`](#arrayexpress_search_experiments) | `sequencing_data` | Searches ArrayExpress functional genomics experiments by keyword, species or platform. |
| [`BioSamples_get_facets`](#BioSamples_get_facets) | `sequencing_data` | Returns BioSamples attribute facets and counts for a text query. |
| [`BioSamples_get_relationships`](#BioSamples_get_relationships) | `sequencing_data` | Retrieves sample-to-sample relationship edges for a BioSample accession. |
| [`BioSamples_get_sample`](#BioSamples_get_sample) | `sequencing_data` | Retrieves a BioSamples record with organism, tissue and linked archive data. |
| [`BioSamples_search`](#BioSamples_search) | `sequencing_data` | Searches EBI BioSamples by text query. |
| [`BioSamples_search_by_filter`](#BioSamples_search_by_filter) | `sequencing_data` | Searches EBI BioSamples by attribute-value filters such as organism or tissue. |
| [`CancerPrognosis_get_study_summary`](#CancerPrognosis_get_study_summary) | `sequencing_data` | Summarizes available molecular profiles, clinical attributes and sample counts of a cancer study. |
| [`CancerPrognosis_search_studies`](#CancerPrognosis_search_studies) | `sequencing_data` | Searches cBioPortal cancer genomics studies by keyword. |
| [`cBioPortal_get_cancer_studies`](#cBioPortal_get_cancer_studies) | `sequencing_data` | Returns a page of the cBioPortal cancer study catalogue. |
| [`cBioPortal_get_gene_panels`](#cBioPortal_get_gene_panels) | `sequencing_data` | Lists gene panels used across cBioPortal studies. |
| [`cBioPortal_get_molecular_profiles`](#cBioPortal_get_molecular_profiles) | `sequencing_data` | Lists molecular profiles such as mutations, CNA and mRNA expression for a study. |
| [`cBioPortal_get_sample_lists`](#cBioPortal_get_sample_lists) | `sequencing_data` | Lists sample lists (cohort definitions) in a cBioPortal study. |
| [`cBioPortal_get_samples`](#cBioPortal_get_samples) | `sequencing_data` | Lists all samples in a cBioPortal study. |
| [`CELLxGENE_download_h5ad`](#CELLxGENE_download_h5ad) | `sequencing_data` | Downloads or returns URIs of original H5AD files from CELLxGENE datasets. |
| [`ChIPAtlas_get_experiment_metadata`](#ChIPAtlas_get_experiment_metadata) | `sequencing_data` | Returns structured metadata for one ChIP-Atlas experiment accession. |
| [`ChIPAtlas_get_experiments`](#ChIPAtlas_get_experiments) | `sequencing_data` | Searches ChIP-Atlas experiment metadata by genome, antigen or cell type. |
| [`ChIPAtlas_search_datasets`](#ChIPAtlas_search_datasets) | `sequencing_data` | Searches ChIP-Atlas for datasets and experiment IDs by antigen or cell type. |
| [`CxGDisc_get_collection`](#CxGDisc_get_collection) | `sequencing_data` | Retrieves a CellxGene Discover collection with its datasets and assay details. |
| [`CxGDisc_list_collections`](#CxGDisc_list_collections) | `sequencing_data` | Lists curated CellxGene Discover single-cell collections. |
| [`CxGDisc_search_datasets`](#CxGDisc_search_datasets) | `sequencing_data` | Searches CellxGene Discover datasets by tissue, disease, organism or cell type. |
| [`DbGaP_get_study`](#DbGaP_get_study) | `sequencing_data` | Gets public metadata for a dbGaP study by phs accession, including conditions and dataset counts. |
| [`DbGaP_search_studies`](#DbGaP_search_studies) | `sequencing_data` | Searches dbGaP studies by title keyword via its FHIR API. |
| [`DDBJ_get_cross_references`](#DDBJ_get_cross_references) | `sequencing_data` | Lists DDBJ records linked to an accession, such as a study's experiments, runs and samples. |
| [`DDBJ_get_entry`](#DDBJ_get_entry) | `sequencing_data` | Retrieves a DDBJ record (BioProject, BioSample, SRA, JGA, GEA, MetaboBank) by accession. |
| [`DDBJ_search_entries`](#DDBJ_search_entries) | `sequencing_data` | Keyword-searches one DDBJ entry type such as BioProject, SRA or GEA, optionally by organism taxid. |
| [`EGA_get_dataset`](#EGA_get_dataset) | `sequencing_data` | Gets public metadata for an EGA dataset accession (platform, sample count, access). |
| [`EGA_get_study`](#EGA_get_study) | `sequencing_data` | Gets public metadata for an EGA study accession. |
| [`EGA_get_study_datasets`](#EGA_get_study_datasets) | `sequencing_data` | Lists EGA dataset accessions belonging to a study. |
| [`ena_get_sequence_xml`](#ena_get_sequence_xml) | `sequencing_data` | Gets XML metadata from ENA for study, sample, run, experiment, analysis or taxon records. |
| [`ENAPortal_count_records`](#ENAPortal_count_records) | `sequencing_data` | Counts ENA records (studies, samples, reads) matching a query. |
| [`ENAPortal_search_runs_fastq`](#ENAPortal_search_runs_fastq) | `sequencing_data` | Searches ENA read runs by taxonomy and library strategy, returning FASTQ download URLs. |
| [`ENAPortal_search_samples`](#ENAPortal_search_samples) | `sequencing_data` | Searches ENA biological samples by taxonomy or text. |
| [`ENAPortal_search_studies`](#ENAPortal_search_studies) | `sequencing_data` | Searches ENA sequencing studies by text or taxonomy. |
| [`ENCODE_get_biosample`](#ENCODE_get_biosample) | `sequencing_data` | Gets ENCODE biosample metadata by accession (type, organism, donor, treatments). |
| [`ENCODE_get_experiment`](#ENCODE_get_experiment) | `sequencing_data` | Gets ENCODE experiment metadata by accession including assay, target, biosample, replicates and files. |
| [`ENCODE_get_file`](#ENCODE_get_file) | `sequencing_data` | Gets ENCODE file metadata by accession including format, output type, assembly and download URL. |
| [`ENCODE_list_files`](#ENCODE_list_files) | `sequencing_data` | Lists ENCODE data files filtered by format, output type or assay. |
| [`ENCODE_search_biosamples`](#ENCODE_search_biosamples) | `sequencing_data` | Searches ENCODE biosamples by organism, type or treatment. |
| [`ENCODE_search_experiments`](#ENCODE_search_experiments) | `sequencing_data` | Searches ENCODE experiments by assay, target, organism or status. |
| [`ENCODE_search_microrna_experiments`](#ENCODE_search_microrna_experiments) | `sequencing_data` | Searches ENCODE microRNA-seq experiments. |
| [`ENCODE_search_rnaseq_experiments`](#ENCODE_search_rnaseq_experiments) | `sequencing_data` | Searches ENCODE RNA-seq experiments by biosample and assay type. |
| [`FourDN_get_download_url`](#FourDN_get_download_url) | `sequencing_data` | Gets the download URL and DRS endpoint for a 4DN data file by accession (requires 4DN access key). |
| [`FourDN_get_experiment_metadata`](#FourDN_get_experiment_metadata) | `sequencing_data` | Gets 4DN experiment metadata: assay type (Hi-C, ATAC-seq, ChIP-seq), biosource, protocol and output files. |
| [`GDC_list_files`](#GDC_list_files) | `sequencing_data` | Lists GDC files filtered by data type (e.g. expression quantification) for analysis pipelines. |
| [`GDC_list_projects`](#GDC_list_projects) | `sequencing_data` | Lists GDC projects (TCGA, TARGET, CGCI) with case and file counts, filterable by program. |
| [`GEO_get_dataset_details`](#GEO_get_dataset_details) | `sequencing_data` | Gets metadata for a GEO series (GSE) including platform, samples and supplementary links. |
| [`geo_get_dataset_info`](#geo_get_dataset_info) | `sequencing_data` | Gets title, summary and metadata for a GEO dataset. |
| [`geo_get_sample_info`](#geo_get_sample_info) | `sequencing_data` | Gets sample characteristics and experimental conditions for a GEO dataset. |
| [`geo_list_supplementary_files`](#geo_list_supplementary_files) | `sequencing_data` | Lists downloadable supplementary/raw files of a GEO GSE or GSM with size, type and download URL. |
| [`geo_search_datasets`](#geo_search_datasets) | `sequencing_data` | Searches GEO datasets by query, organism, study type or platform and returns IDs and basic info. |
| [`GEO_search_rnaseq_datasets`](#GEO_search_rnaseq_datasets) | `sequencing_data` | Searches NCBI GEO for RNA-seq datasets by keyword, tissue or gene. |
| [`GSA_get_accession`](#GSA_get_accession) | `sequencing_data` | Looks up a China NGDC Genome Sequence Archive CRA accession returning title and linked BioProject/BioSample. |
| [`GTEx_get_sample_info`](#GTEx_get_sample_info) | `sequencing_data` | Gets GTEx sample and subject metadata: tissue, donor age/sex, ischemic time, RIN and QC metrics. |
| [`hca_get_file_manifest`](#hca_get_file_manifest) | `sequencing_data` | Retrieves downloadable file list for a Human Cell Atlas project by project ID. |
| [`hca_search_projects`](#hca_search_projects) | `sequencing_data` | Searches Human Cell Atlas DCP single-cell projects by organ and disease. |
| [`IGSR_list_data_collections`](#IGSR_list_data_collections) | `sequencing_data` | Lists IGSR data collections (1000 Genomes phases, 30x resequencing, HGSVC, HPRC). |
| [`IGSR_search_samples`](#IGSR_search_samples) | `sequencing_data` | Searches 1000 Genomes samples by population, data collection or name with sex and BioSample IDs. |
| [`IntOGen_list_cohorts`](#IntOGen_list_cohorts) | `sequencing_data` | Lists IntOGen cancer cohorts with data source, sample counts, tumor type and age group, optionally filtered by cancer type. |
| [`MGnify_get_samples`](#MGnify_get_samples) | `sequencing_data` | Lists MGnify samples with BioSample ID, location, host and environment metadata. |
| [`MGnify_get_study_detail`](#MGnify_get_study_detail) | `sequencing_data` | Returns MGnify study details with abstract, BioProject, sample count and biomes. |
| [`MGnify_get_taxonomy`](#MGnify_get_taxonomy) | `sequencing_data` | Returns taxonomic composition and abundances from a MGnify analysis using SSU or LSU rRNA markers. |
| [`MGnify_list_analyses`](#MGnify_list_analyses) | `sequencing_data` | Lists MGnify analyses for a study accession with taxonomic and functional outputs. |
| [`MGnify_list_analysis_downloads`](#MGnify_list_analysis_downloads) | `sequencing_data` | Lists MGnify analysis output files such as predicted CDS and taxonomy tables with download URLs. |
| [`MGnify_search_studies`](#MGnify_search_studies) | `sequencing_data` | Searches MGnify metagenomics/microbiome studies by biome or keyword. |
| [`MPD_get_phenotype_data`](#MPD_get_phenotype_data) | `sequencing_data` | Searches ENCODE experiment records that mention a given mouse strain. |
| [`NCBI_SRA_get_download_urls`](#NCBI_SRA_get_download_urls) | `sequencing_data` | Returns FTP, AWS S3 and web URLs for downloading SRA run data. |
| [`NCBI_SRA_get_run_info`](#NCBI_SRA_get_run_info) | `sequencing_data` | Returns SRA run metadata including platform, instrument, library strategy, layout and sample details. |
| [`NCBI_SRA_link_to_biosample`](#NCBI_SRA_link_to_biosample) | `sequencing_data` | Links SRA runs to BioSample records with sample attributes, tissue, cell line and donor information. |
| [`NCBI_SRA_locate_run_files`](#NCBI_SRA_locate_run_files) | `sequencing_data` | Resolves cloud download locations with file size and md5 for an SRA run via the SRA Data Locator. |
| [`NCBI_SRA_search_runs`](#NCBI_SRA_search_runs) | `sequencing_data` | Searches NCBI SRA runs by study, organism, strategy, platform or library source. |
| [`Nextstrain_get_dataset`](#Nextstrain_get_dataset) | `sequencing_data` | Returns metadata and phylogenetic summary for a Nextstrain pathogen dataset, including sequence counts. |
| [`Nextstrain_list_datasets`](#Nextstrain_list_datasets) | `sequencing_data` | Lists Nextstrain pathogen phylogenetic datasets grouped by pathogen. |
| [`OSDR_get_study`](#OSDR_get_study) | `sequencing_data` | Returns metadata for a NASA OSDR space biology study including organism, mission, factors and assay type. |
| [`OSDR_list_files`](#OSDR_list_files) | `sequencing_data` | Lists downloadable raw and processed data files for a NASA OSDR study. |
| [`OSDR_search_studies`](#OSDR_search_studies) | `sequencing_data` | Searches NASA OSDR space biology studies across transcriptomics, proteomics and other assays. |
| [`Pathoplexus_get_sequence_details`](#Pathoplexus_get_sequence_details) | `sequencing_data` | Lists per-sequence metadata records (accession, country, date, host, lineage) for a Pathoplexus organism. |
| [`PharmacoDB_get_molecular_profiling`](#PharmacoDB_get_molecular_profiling) | `sequencing_data` | Lists which molecular profiling data types (RNA-seq, mutation, CNV) exist for a cell line and in which dataset. |
| [`Progenetix_list_cohorts`](#Progenetix_list_cohorts) | `sequencing_data` | Lists Progenetix cancer CNV cohorts (e.g. TCGA) with IDs, sizes and data types. |
| [`Progenetix_search_biosamples`](#Progenetix_search_biosamples) | `sequencing_data` | Searches Progenetix tumor biosamples by NCIt disease code, returning IDs, histology and stage. |
| [`SCP_get_study`](#SCP_get_study) | `sequencing_data` | Retrieves one Broad Single Cell Portal study by accession with its full description and metadata. |
| [`SCP_list_studies`](#SCP_list_studies) | `sequencing_data` | Lists public Single Cell Portal studies sorted by cell count, optionally filtered by minimum cells or keyword. |
| [`SCP_search_studies`](#SCP_search_studies) | `sequencing_data` | Searches Broad Single Cell Portal studies by keyword and returns accession, title, cell and gene counts and URL. |
| [`SCXA_list_experiments`](#SCXA_list_experiments) | `sequencing_data` | Lists curated single-cell RNA-seq experiments in the EBI Single Cell Expression Atlas with species, cell counts, technology and experimental factors. |
| [`SRA_get_experiment`](#SRA_get_experiment) | `sequencing_data` | Gets metadata for an SRA sequencing experiment by accession, including organism, platform, library details, study and run statistics. |
| [`SRA_search_experiments`](#SRA_search_experiments) | `sequencing_data` | Searches NCBI SRA for sequencing experiments by keyword, organism, library strategy or platform. |
| [`UCSCCellBrowser_get_dataset`](#UCSCCellBrowser_get_dataset) | `sequencing_data` | Gets details of a UCSC Cell Browser single-cell dataset including abstract, organisms, body parts, diseases and child datasets. |
| [`UCSCCellBrowser_list_facets`](#UCSCCellBrowser_list_facets) | `sequencing_data` | Lists distinct organism, body part or disease facet values with dataset counts in the UCSC Cell Browser catalog. |
| [`UCSCCellBrowser_search_datasets`](#UCSCCellBrowser_search_datasets) | `sequencing_data` | Searches the UCSC Cell Browser catalog of single-cell datasets by organism, body part, disease or keyword. |

## `arrayexpress_get_experiment` {#arrayexpress_get_experiment}

Retrieves ArrayExpress experiment metadata by accession (E-MTAB, E-GEOD).

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `experiment_id`
- **Keywords:** arrayexpress, functional genomics, experiment, e-mtab
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `experiment_id` | string | yes | ArrayExpress experiment ID (e.g., 'E-MTAB-1234') |

??? note "Upstream description"

    Get detailed information about a specific ArrayExpress experiment by its accession ID (E-GEOD, E-MTAB, etc.). Returns complete experiment metadata including title, description, protocols, experimental design, and publication information from the original ArrayExpress database.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'arrayexpress_get_experiment',
     'arguments': {'experiment_id': '...'}},
)
```

## `arrayexpress_get_experiment_files` {#arrayexpress_get_experiment_files}

Lists raw, processed and MAGE-TAB files with download paths for an ArrayExpress experiment.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `experiment_id`
- **Keywords:** raw data, processed data, mage-tab, download
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `experiment_id` | string | yes | ArrayExpress experiment ID |

??? note "Upstream description"

    Get list of files associated with an ArrayExpress experiment. Returns file metadata including file names, types, sizes, and download paths for raw data files, processed data files, and MAGE-TAB annotation files.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'arrayexpress_get_experiment_files',
     'arguments': {'experiment_id': '...'}},
)
```

## `arrayexpress_get_experiment_samples` {#arrayexpress_get_experiment_samples}

Returns sample characteristics and experimental factors for an ArrayExpress experiment.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `experiment_id`
- **Keywords:** samples, sdrf, experimental factors, mage-tab
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `experiment_id` | string | yes | ArrayExpress experiment ID |

??? note "Upstream description"

    Get sample information for an ArrayExpress experiment. Returns sample metadata including sample names, characteristics (e.g., organism part, developmental stage), experimental factors, and technical details from the original ArrayExpress MAGE-TAB annotations.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'arrayexpress_get_experiment_samples',
     'arguments': {'experiment_id': '...'}},
)
```

## `arrayexpress_search_experiments` {#arrayexpress_search_experiments}

Searches ArrayExpress functional genomics experiments by keyword, species or platform.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** experiment search, microarray, rna-seq, species
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `keywords` | string | no | Search keywords (e.g., disease name, tissue type, experimental condition) |
| `species` | string | no | Species name (e.g., 'Homo sapiens', 'Mus musculus') |
| `array` | string | no | Array platform name |
| `limit` | integer | no | Maximum number of results (default: 10, max: 100) Default `10` |
| `offset` | integer | no | Offset for pagination (default: 0) Default `0` |

??? note "Upstream description"

    Search ArrayExpress experiments - the original functional genomics database from EMBL-EBI. Searches microarray and RNA-seq experiments by keywords, species, or array platform. Returns experiment metadata including accession IDs (E-GEOD, E-MTAB), titles, descriptions, and study information. Note: Some newer ArrayExpress data is now hosted in BioStudies - use BioStudies tools if an accession is not found here.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'arrayexpress_search_experiments', 'arguments': {}},
)
```

## `BioSamples_get_facets` {#BioSamples_get_facets}

Returns BioSamples attribute facets and counts for a text query.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `text`
- **Keywords:** biosamples, facets, sample metadata, sra accession
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `text` | string | yes | Free-text query to compute facets over. Examples: 'cancer', 'liver', 'Homo sapiens'. |
| `max_values` | integer \| null | no | Maximum number of top values to return per facet (1-50, default 10). |

??? note "Upstream description"

    Faceted aggregation discovery over EBI BioSamples: given a free-text query, return the available attribute facets (organism, status, SRA accession, external reference, release date range, ...) with sample counts plus the top values per facet, so a scientist learns which filterable attribute keys/values exist before running a structured filter search (BioSamples_search_by_filter). Also returns external-reference data facets (e.g. ENA). Example: text='cancer' returns organism (count 59585), SRA accession (54248), status (59590), external reference (54524), release date range (59590), and ENA external-reference facet (53344). No authentication required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BioSamples_get_facets', 'arguments': {'text': '...'}},
)
```

## `BioSamples_get_relationships` {#BioSamples_get_relationships}

Retrieves sample-to-sample relationship edges for a BioSample accession.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** sample lineage, derived from, provenance
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | BioSamples accession. Formats: SAMEA* (EBI), SAMN* (NCBI), SAMD* (DDBJ), SAMEG* (sample group). Example: 'SAMEA4451312'. |

??? note "Upstream description"

    Retrieve the relationship/provenance graph of a BioSample: the sample-to-sample lineage edges (e.g. 'derived from', 'has member', 'child of', 'same as', 'recurated into') that the standard BioSamples_get_sample parser drops. Each edge reports source, type, target, and the direction relative to the queried sample, so processed-sample -&gt; source-tissue and sample-group -&gt; member-sample provenance becomes reachable. Example: SAMEA4451312 returns edges \[SAMEA4451312 derived from SAMEA4451117\] and \[SAMEG315830 has member SAMEA4451312\]. No authentication required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BioSamples_get_relationships',
     'arguments': {'accession': '...'}},
)
```

## `BioSamples_get_sample` {#BioSamples_get_sample}

Retrieves a BioSamples record with organism, tissue and linked archive data.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** biosample accession, samea, ena link, sample metadata
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | BioSamples accession. Formats: SAMEA* (EBI), SAMN* (NCBI), SAMD* (DDBJ). Examples: 'SAMEA104228123', 'SAMD00000603'. |

??? note "Upstream description"

    Get detailed metadata for a specific biological sample from the EBI BioSamples database by accession. Returns sample name, organism, tissue, disease, experimental context, and links to associated data archives (ENA, ArrayExpress, EVA). BioSamples contains 60+ million samples from all major sequence and functional genomics archives. Example: SAMEA104228123 returns a metagenomic sample from the Qiita microbiome project.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BioSamples_get_sample',
     'arguments': {'accession': '...'}},
)
```

## `BioSamples_search` {#BioSamples_search}

Searches EBI BioSamples by text query.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `query`
- **Keywords:** sample search, biosamples, organism
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Free-text search query. Examples: 'breast cancer', 'liver tissue', 'CRISPR screen', 'Mus musculus brain'. |
| `limit` | integer \| null | no | Maximum results to return (1-50, default 10). |
| `page` | integer \| null | no | Page number (0-indexed, default 0). |

??? note "Upstream description"

    Search the EBI BioSamples database by text query across 60+ million biological samples. Returns matching samples with accessions, names, organisms, and key characteristics. Useful for finding source material metadata for sequencing, proteomics, and other experiments. Example: searching 'breast cancer' returns thousands of tumor samples with tissue/disease annotations.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BioSamples_search', 'arguments': {'query': '...'}},
)
```

## `BioSamples_search_by_filter` {#BioSamples_search_by_filter}

Searches EBI BioSamples by attribute-value filters such as organism or tissue.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `attribute`, `value`
- **Keywords:** attribute filter, sample metadata, tissue
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `attribute` | string | yes | Attribute name to filter on. Common attributes: 'organism', 'tissue', 'disease', 'cell type', 'sex', 'age'. |
| `value` | string | yes | Value to match for the attribute. Examples: 'Homo sapiens', 'liver', 'melanoma', 'female'. |
| `limit` | integer \| null | no | Maximum results to return (1-50, default 10). |

??? note "Upstream description"

    Search the EBI BioSamples database with structured attribute filters. More precise than text search - filters by specific metadata fields like organism, tissue, or disease. Example: attribute='organism', value='Homo sapiens' returns human samples; attribute='disease', value='melanoma' returns melanoma-specific samples.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'BioSamples_search_by_filter',
     'arguments': {'attribute': '...', 'value': '...'}},
)
```

## `CancerPrognosis_get_study_summary` {#CancerPrognosis_get_study_summary}

Summarizes available molecular profiles, clinical attributes and sample counts of a cancer study.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** cancer study, molecular profiles, tcga, cnv
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `get_study_summary` Default `get_study_summary` |
| `cancer` | string | no | TCGA cancer type abbreviation (e.g., 'BRCA', 'LUAD', 'COADREAD') or cBioPortal study ID. Use cancer_type as alias. |
| `cancer_type` | string | no | Alias for cancer. TCGA cancer type abbreviation (e.g., 'BRCA', 'LUAD', 'COADREAD'). |
| `study_id` | string | no | Alias for cancer. cBioPortal study ID (e.g., 'brca_tcga', 'luad_tcga'). |

??? note "Upstream description"

    Get summary information for a cancer study including available molecular profiles (mutations, expression, CNV), survival-related clinical attributes, and sample counts. Useful for determining what data is available before querying expression or survival data. Lists all 33 supported TCGA cancer type abbreviations.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CancerPrognosis_get_study_summary', 'arguments': {}},
)
```

## `CancerPrognosis_search_studies` {#CancerPrognosis_search_studies}

Searches cBioPortal cancer genomics studies by keyword.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** cancer genomics, tcga, study search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `search_studies` Default `search_studies` |
| `keyword` | string | no | Search keyword (e.g., 'breast', 'lung', 'TCGA', 'melanoma', 'glioblastoma') |
| `limit` | integer \| null | no | Maximum number of results to return (default 20, max 100) |
| `cancer_type` | string | no | Cancer type to search for. Alias for keyword. |
| `cancer` | string | no | Cancer type to search for. Alias for keyword. |
| `query` | string | no | Search query. Alias for keyword (e.g., 'lung', 'breast', 'TCGA'). |

??? note "Upstream description"

    Search cBioPortal for cancer genomics studies by keyword. Find studies by cancer type, institution, or data type. Returns study IDs that can be used with other CancerPrognosis tools. Covers 400+ studies including all TCGA projects, AACR GENIE, MSK-IMPACT, and institutional cohorts.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CancerPrognosis_search_studies', 'arguments': {}},
)
```

## `cBioPortal_get_cancer_studies` {#cBioPortal_get_cancer_studies}

Returns a page of the cBioPortal cancer study catalogue.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** cbioportal, cancer studies, tcga
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `limit` | integer | no | Number of studies to return in this page. The default of 20 covers only a small fraction of the catalogue; set it to the `total_available` reported in the response to retrieve every study. Default `20` |
| `offset` | integer | no | 0-based index of the first study to return, for paging through the catalogue. Applied client-side because the cBioPortal /studies endpoint ignores its own pageNumber parameter. Default `0` |

??? note "Upstream description"

    Get one page of the cBioPortal study catalogue. This returns a SLICE, not the whole catalogue: `limit` defaults to 20 while cBioPortal currently holds several hundred studies. The response always reports `total_available` (the full catalogue size) alongside `count` (studies actually returned), and sets `truncated: true` with a `truncation_note` whenever studies were left behind. Never conclude a study is missing from cBioPortal from a truncated page — raise `limit` to `total_available`, or page with `offset`, before drawing that conclusion.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'cBioPortal_get_cancer_studies', 'arguments': {}},
)
```

## `cBioPortal_get_gene_panels` {#cBioPortal_get_gene_panels}

Lists gene panels used across cBioPortal studies.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** gene panels, targeted sequencing, coverage
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `page_size` | integer | no | Number of panels to return Default `50` |

??? note "Upstream description"

    Get gene panels used in cBioPortal studies. Gene panels define which genes were sequenced in a study. Use to understand coverage before querying mutation data. Returns one page: `page_size` defaults to 50, and the response reports `total_available` (all panels matching) next to `count` (panels returned), plus `truncated`/`truncation_note` when panels were left behind.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'cBioPortal_get_gene_panels', 'arguments': {}},
)
```

## `cBioPortal_get_molecular_profiles` {#cBioPortal_get_molecular_profiles}

Lists molecular profiles such as mutations, CNA and mRNA expression for a study.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `study_id`
- **Keywords:** molecular profiles, mutation profile, mrna, cna
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `study_id` | string | yes | Cancer study ID |

??? note "Upstream description"

    Get molecular profiles for a cancer study. Molecular profiles include mutation data, copy number alterations, mRNA expression, etc.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'cBioPortal_get_molecular_profiles',
     'arguments': {'study_id': '...'}},
)
```

## `cBioPortal_get_sample_lists` {#cBioPortal_get_sample_lists}

Lists sample lists (cohort definitions) in a cBioPortal study.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `study_id`
- **Keywords:** sample lists, cohort, sequenced samples
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `study_id` | string | yes | Cancer study ID |

??? note "Upstream description"

    Get all sample lists (cohort definitions) in a cancer study. Sample lists define groups like 'samples with mutations', 'samples with CNA data'. Use for targeted queries.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'cBioPortal_get_sample_lists',
     'arguments': {'study_id': '...'}},
)
```

## `cBioPortal_get_samples` {#cBioPortal_get_samples}

Lists all samples in a cBioPortal study.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `study_id`
- **Keywords:** samples, tumor samples, cohort size
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `study_id` | string | yes | Cancer study ID |
| `page_size` | integer | no | Maximum number of samples to return. Defaults high enough to return every sample in a study; a smaller value silently truncates the cohort, which will corrupt any frequency denominator computed from it. Default `100000` |
| `page_number` | integer | no | 0-based page index, used with `page_size` to page through very large studies. Default `0` |

??? note "Upstream description"

    Get all samples in a cancer study. Returns the full cohort by default; the response is a bare array with no truncation marker, so lower `page_size` only when you deliberately want a sample. Use this count as the denominator for mutation-frequency calculations.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'cBioPortal_get_samples',
     'arguments': {'study_id': '...'}},
)
```

## `CELLxGENE_download_h5ad` {#CELLxGENE_download_h5ad}

Downloads or returns URIs of original H5AD files from CELLxGENE datasets.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `dataset_id`
- **Keywords:** h5ad, anndata, single-cell dataset, download
- **Needs API keys:** `CELLXGENE_CENSUS_PACKAGE_INSTALLED`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `download_h5ad` Default `download_h5ad` |
| `dataset_id` | string | yes | CELLxGENE dataset identifier (required) |
| `output_path` | string | no | Local path to save H5AD file (omit to get URI only) |
| `census_version` | string | no | Census version to query. 'stable' (recommended, Long-Term Support release), 'latest' (newest data, may change), or specific date 'YYYY-MM-DD' (for reproducibility). Default 'stable' is best for production analyses. Default `stable` |

??? note "Upstream description"

    Download original H5AD (HDF5-based AnnData) files from CELLxGENE datasets or get their URIs. Access unprocessed data from individual studies with complete metadata and analysis results. Requires dataset_id from CELLxGENE Discover. Note: Files can be large (GBs). Use for: accessing raw data, reproducibility, offline analysis, custom processing pipelines. Note: requires the cellxgene-census Python package; install with `pip install cellxgene-census` for live access.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CELLxGENE_download_h5ad',
     'arguments': {'dataset_id': '...'}},
)
```

## `ChIPAtlas_get_experiment_metadata` {#ChIPAtlas_get_experiment_metadata}

Returns structured metadata for one ChIP-Atlas experiment accession.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `experiment_id`
- **Keywords:** srx, chip-seq experiment, antigen, cell type
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | One of `get_experiment_metadata` Default `get_experiment_metadata` |
| `experiment_id` | string | yes | Experiment accession in SRX/DRX/ERX format (e.g. 'SRX080331'). |
| `genome` | string | no | Optional: restrict to one genome assembly (an experiment can have records for several, e.g. hg38 and hg19). One of `hg38`, `hg19`, `mm10`, `mm9`, `rn6`, `dm6`, `dm3`, `ce11`, `ce10`, `sacCer3` |

??? note "Upstream description"

    ChIP-Atlas single-experiment structured metadata: given one experiment accession (SRX/DRX/ERX), return its antigen class and antigen, cell-type class and cell type, title, and the full parsed sample attributes (antibody, source, cell line, etc.). One accession can map to multiple genome assemblies; filter with the optional 'genome' parameter. Use for: resolving what an SRX experiment actually measured, looking up the antibody/cell line behind a peak set, validating experiment IDs returned by other ChIP-Atlas tools.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ChIPAtlas_get_experiment_metadata',
     'arguments': {'experiment_id': '...'}},
)
```

## `ChIPAtlas_get_experiments` {#ChIPAtlas_get_experiments}

Searches ChIP-Atlas experiment metadata by genome, antigen or cell type.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** chip-seq, experiments, srx, antigen
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | One of `get_experiment_list` Default `get_experiment_list` |
| `genome` | string | no | Filter by genome assembly One of `hg38`, `hg19`, `mm10`, `mm9`, `rn6`, `dm6`, `dm3`, `ce11`, `ce10`, `sacCer3` |
| `antigen` | string | no | Filter by antigen/protein name (e.g., 'CTCF', 'H3K4me3') |
| `cell_type` | string | no | Filter by cell type (e.g., 'K-562', 'HepG2', 'CD4') |
| `limit` | integer | no | Maximum results to return Default `100` |

??? note "Upstream description"

    Search ChIP-Atlas experiment metadata including experiment IDs, antigens, cell types, and processing statistics. Filter by genome, antigen (TF/histone), or cell type. Returns experiment details with SRX/ERX/DRX IDs. Use for: finding relevant experiments, getting experiment metadata, dataset discovery, batch analysis planning.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ChIPAtlas_get_experiments', 'arguments': {}},
)
```

## `ChIPAtlas_search_datasets` {#ChIPAtlas_search_datasets}

Searches ChIP-Atlas for datasets and experiment IDs by antigen or cell type.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** chip-atlas, dataset search, antigen list, cell type list
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | One of `search_datasets` Default `search_datasets` |
| `antigen` | string | no | Search by antigen/protein name (e.g., 'CTCF', 'H3K27ac') |
| `cell_type` | string | no | Search by cell type (e.g., 'K-562', 'HeLa', 'MCF-7') |
| `genome` | string | no | Genome assembly One of `hg38`, `hg19`, `mm10`, `mm9`, `rn6`, `dm6`, `dm3`, `ce11`, `ce10`, `sacCer3` Default `hg38` |

??? note "Upstream description"

    Search ChIP-Atlas by antigen or cell type to find available datasets. Returns number of experiments and experiment IDs. Use antigenList for finding TF/histone datasets across cell types. Use celltypeList for finding all factors in specific cell type. Use for: dataset discovery, planning comparative analysis, checking data availability.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ChIPAtlas_search_datasets', 'arguments': {}},
)
```

## `CxGDisc_get_collection` {#CxGDisc_get_collection}

Retrieves a CellxGene Discover collection with its datasets and assay details.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `collection_id`
- **Keywords:** collection, single-cell datasets, cellxgene
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `collection_id` | string | yes | Collection UUID from CellxGene Discovery (e.g., '48259aa8-f168-4bf5-b797-af8e88da6637' for Human Breast Cell Atlas). |

??? note "Upstream description"

    Get detailed information about a specific CellxGene Discovery collection, including all datasets with tissue, disease, organism, cell type, and assay details. Example: collection '48259aa8-f168-4bf5-b797-af8e88da6637' is the Human Breast Cell Atlas with 5 datasets (stroma, epithelial, integrated) totaling 2.8M cells. Each dataset includes cell count, tissue labels, disease status, and available assay types.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CxGDisc_get_collection',
     'arguments': {'collection_id': '...'}},
)
```

## `CxGDisc_list_collections` {#CxGDisc_list_collections}

Lists curated CellxGene Discover single-cell collections.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** collections, scrna-seq, cellxgene
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `limit` | integer | no | Maximum number of collections to return (1-100). Default: 20. |

??? note "Upstream description"

    List curated single-cell RNA-seq collections from CZI CellxGene Discovery. Each collection is a published dataset group with DOI, typically corresponding to a paper. The platform hosts 350+ public collections spanning diverse tissues, diseases, and organisms with cell-type annotations. Returns collection name, DOI, dataset count, total cell count, and curator info.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CxGDisc_list_collections', 'arguments': {}},
)
```

## `CxGDisc_search_datasets` {#CxGDisc_search_datasets}

Searches CellxGene Discover datasets by tissue, disease, organism or cell type.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** dataset search, scrna-seq, tissue, disease
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `tissue` | string | no | Filter by tissue name (e.g., 'lung', 'brain', 'heart', 'kidney', 'liver'). Case-insensitive partial match. |
| `disease` | string | no | Filter by disease (e.g., 'cancer', 'COVID', 'Alzheimer', 'normal'). Case-insensitive partial match. |
| `organism` | string | no | Filter by organism (e.g., 'Homo sapiens', 'Mus musculus'). Case-insensitive partial match. |
| `cell_type` | string | no | Filter by cell type (e.g., 'T cell', 'neuron', 'macrophage', 'epithelial'). Case-insensitive partial match. |
| `limit` | integer | no | Maximum datasets to return (1-100). Default: 20. |

??? note "Upstream description"

    Search CellxGene Discovery single-cell RNA-seq datasets by tissue, disease, organism, or cell type. The platform hosts 2,000+ datasets with detailed cell-type annotations. Supports flexible text-based filtering: 'lung' returns 134 lung datasets, 'brain' returns 144, 'T cell' returns 667, 'cancer' returns 193. Results sorted by cell count (largest first).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'CxGDisc_search_datasets', 'arguments': {}},
)
```

## `DbGaP_get_study` {#DbGaP_get_study}

Gets public metadata for a dbGaP study by phs accession, including conditions and dataset counts.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `phs_id`
- **Keywords:** dbgap, phs accession, genotype-phenotype, controlled access, study
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `phs_id` | string | yes | dbGaP study accession, e.g. 'phs000681'. |

??? note "Upstream description"

    Retrieve one dbGaP study's full public metadata by its phs accession: description, medical conditions studied, release date, and phenotype/molecular dataset and variable counts. Example: phs_id='phs000681' returns the Diabetes Prevention Program Genetics study with 4 phenotype datasets and 15 variables. Metadata only -- the underlying genotype/phenotype data requires Data Access Committee approval. Use DbGaP_search_studies to find a phs_id.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DbGaP_get_study', 'arguments': {'phs_id': '...'}},
)
```

## `DbGaP_search_studies` {#DbGaP_search_studies}

Searches dbGaP studies by title keyword via its FHIR API.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `query`
- **Keywords:** dbgap, genotype data, controlled access, study search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Study title keyword, e.g. 'diabetes', 'autism'. |
| `limit` | integer \| null | no | Max studies to return, 1-100. Default 25. |

??? note "Upstream description"

    Search dbGaP (NCBI's controlled-access US genotype-phenotype archive, the American counterpart to EGA) for studies by title keyword, via its FHIR API. Example: query='diabetes' returns the Diabetes Prevention Program Genetics study. Returns each study's phs_id for use with DbGaP_get_study.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DbGaP_search_studies', 'arguments': {'query': '...'}},
)
```

## `DDBJ_get_cross_references` {#DDBJ_get_cross_references}

Lists DDBJ records linked to an accession, such as a study's experiments, runs and samples.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** ddbj, sra, bioproject, biosample, cross-reference
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | DDBJ accession, e.g. 'DRP000001'. |
| `entry_type` | string \| null | no | Override the inferred entry type. One of: bioproject, biosample, sra-study, sra-experiment, sra-run, sra-sample, jga-study. |
| `reference_type` | string \| null | no | Only return links of this type, e.g. 'sra-run', 'bioproject'. |
| `limit` | integer \| null | no | Maximum links to return (default 25, max 200). |

??? note "Upstream description"

    List the DDBJ records linked to an accession, so you can walk a study to its experiments, runs, and samples, or back to its BioProject. Returns linked accessions with their type and URL, plus a count of links per type. Filter with reference_type, e.g. 'sra-run'. Studies can link thousands of runs, so results are capped by limit.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DDBJ_get_cross_references',
     'arguments': {'accession': '...'}},
)
```

## `DDBJ_get_entry` {#DDBJ_get_entry}

Retrieves a DDBJ record (BioProject, BioSample, SRA, JGA, GEA, MetaboBank) by accession.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** ddbj, insdc, sra, jga, accession
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | DDBJ accession, e.g. 'DRP000001', 'PRJDB3490', 'JGAS000001', 'E-GEAD-1000', 'MTBKS102'. |
| `entry_type` | string \| null | no | Override the inferred entry type. One of: bioproject, biosample, sra-study, sra-experiment, sra-run, sra-sample, jga-study, gea, metabobank. One of `None`, `bioproject`, `biosample`, `sra-study`, `sra-experiment`, `sra-run`, `sra-sample`, `jga-study`, `gea`, `metabobank` |

??? note "Upstream description"

    Retrieve a DNA Data Bank of Japan record by accession. DDBJ is the third INSDC member alongside NCBI/GenBank and EMBL-EBI/ENA, and also hosts JGA (Japan's controlled-access human archive), GEA (expression), and MetaboBank (metabolomics). The entry type is inferred from the accession prefix (PRJD*=bioproject, SAMD*=biosample, DRP*=sra-study, DRX*=sra-experiment, DRR*=sra-run, DRS*=sra-sample, JGAS*=jga-study, E-GEAD*=gea, MTBKS*=metabobank) unless you pass entry_type. Example: accession='DRP000001'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DDBJ_get_entry', 'arguments': {'accession': '...'}},
)
```

## `DDBJ_search_entries` {#DDBJ_search_entries}

Keyword-searches one DDBJ entry type such as BioProject, SRA or GEA, optionally by organism taxid.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `entry_type`, `keywords`
- **Keywords:** ddbj, sra, gea, bioproject, search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `entry_type` | string | yes | Which DDBJ record type to search. One of `bioproject`, `biosample`, `sra-study`, `sra-experiment`, `sra-run`, `sra-sample`, `jga-study`, `gea`, `metabobank` |
| `keywords` | string | yes | Free-text search terms, e.g. 'daptomycin resistance'. |
| `organism_taxid` | string \| integer \| null | no | Optional NCBI taxonomy ID filter, e.g. 9606 (human). |
| `limit` | integer \| null | no | Max entries to return, 1-100. Default 25. |

??? note "Upstream description"

    Keyword search one DDBJ entry type (bioproject, biosample, sra-study, sra-experiment, sra-run, sra-sample, jga-study, gea, metabobank), optionally filtered by organism NCBI taxonomy ID. Example: entry_type='gea', keywords='daptomycin resistance' finds gene-expression studies on antibiotic resistance. Returns each match's identifier for use with DDBJ_get_entry or DDBJ_get_cross_references.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'DDBJ_search_entries',
     'arguments': {'entry_type': 'bioproject', 'keywords': '...'}},
)
```

## `EGA_get_dataset` {#EGA_get_dataset}

Gets public metadata for an EGA dataset accession (platform, sample count, access).

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** ega, controlled access, genotyping, dataset
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | EGA dataset accession, e.g. 'EGAD00000000001'. |

??? note "Upstream description"

    Resolve an EGA dataset accession to its public metadata: technology platform, sample count, and access policy. Example: accession='EGAD00000000001' returns an Affymetrix 500K genotyping dataset of 1,504 samples under controlled access. Metadata only -- the underlying data requires Data Access Committee approval via the returned policy_accession_id.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EGA_get_dataset', 'arguments': {'accession': '...'}},
)
```

## `EGA_get_study` {#EGA_get_study}

Gets public metadata for an EGA study accession.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** ega, study, genome-phenome, controlled access
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | EGA study accession, e.g. 'EGAS00000000001'. |

??? note "Upstream description"

    Resolve an EGA (European Genome-phenome Archive) study accession to its public metadata: title, description, study type, and publication. Example: accession='EGAS00000000001' returns the WTCCC bipolar disorder case-control study. Metadata only -- the underlying controlled-access sequence/genotype data requires separate Data Access Committee approval. Use EGA_get_study_datasets to find the study's dataset accessions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EGA_get_study', 'arguments': {'accession': '...'}},
)
```

## `EGA_get_study_datasets` {#EGA_get_study_datasets}

Lists EGA dataset accessions belonging to a study.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** ega, datasets, study, egad
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | EGA study accession, e.g. 'EGAS00000000001'. |

??? note "Upstream description"

    List the dataset accessions belonging to one EGA study. Example: accession='EGAS00000000001' returns its constituent EGAD dataset(s) with technology and sample count. Use EGA_get_study to find a study accession first, then EGA_get_dataset for full details on any returned dataset.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'EGA_get_study_datasets',
     'arguments': {'accession': '...'}},
)
```

## `ena_get_sequence_xml` {#ena_get_sequence_xml}

Gets XML metadata from ENA for study, sample, run, experiment, analysis or taxon records.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** ena, xml, study, run, experiment
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | ENA metadata record accession. Supported types: Study (ERP*, SRP*, PRJ*), Sample (ERS*, SRS*, SAM*), Run (ERR*, SRR*), Experiment (ERX*, SRX*), Analysis (ERZ*, SRZ*), Taxon (numeric ID). NOT for sequence accessions (U\#\#\#\#\#, M\#\#\#\#\#, etc.). Examples: 'SAMEA2620084' (sample), 'ERX2484999' (experiment), '562' (taxon) |

??? note "Upstream description"

    Get metadata in XML format from ENA for Study, Sample, Run, Experiment, Analysis, or Taxon records. XML format is NOT available for sequence records - use ena_get_sequence_fasta or ena_get_sequence_embl for sequence accessions. Supported accession types: Study (ERP*, SRP*, PRJ*), Sample (ERS*, SRS*, SAM*), Run (ERR*, SRR*), Experiment (ERX*, SRX*), Analysis (ERZ*, SRZ*), Taxon IDs. Returns XML-formatted metadata record with submitter-provided information, cross-references, and FTP locations.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ena_get_sequence_xml', 'arguments': {'accession': '...'}},
)
```

## `ENAPortal_count_records` {#ENAPortal_count_records}

Counts ENA records (studies, samples, reads) matching a query.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `query`
- **Keywords:** ena, record count, studies, reads
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | ENA search query. Examples: 'description="cancer"', 'tax_tree(9606)', 'description="CRISPR" AND tax_tree(10090)'. |
| `result_type` | string \| null | no | Type of records to count. Options: 'study', 'sample', 'read_run', 'analysis', 'assembly'. Default: 'study'. |

??? note "Upstream description"

    Count the number of records matching a query in the European Nucleotide Archive (ENA). Supports counting studies, samples, reads, and other result types. Useful for estimating dataset sizes before downloading. Example: counting studies with 'description="transcriptomics"' returns ~9,000+, showing the breadth of transcriptomics data in ENA.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENAPortal_count_records', 'arguments': {'query': '...'}},
)
```

## `ENAPortal_search_runs_fastq` {#ENAPortal_search_runs_fastq}

Searches ENA read runs by taxonomy and library strategy, returning FASTQ download URLs.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `query`
- **Keywords:** ena, fastq, read run, library strategy
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | ENA read_run query. Examples: 'tax_eq(1280)' (S. aureus, exact), 'tax_tree(590)' (Salmonella incl. descendants), 'tax_eq(1280) AND library_strategy="WGS"'. Plain keywords are auto-wrapped as description="...". |
| `limit` | integer \| null | no | Maximum runs to return (1-100, default 10). |
| `fields` | string \| null | no | Comma-separated ENA run fields. Default: 'run_accession,fastq_ftp,submitted_ftp,sra_ftp,fastq_md5,fastq_bytes,library_strategy,instrument_platform,read_count,scientific_name'. |

??? note "Upstream description"

    Search the European Nucleotide Archive (ENA) for sequencing RUNS (read_run) by taxonomy and library strategy, returning ready-to-use FASTQ download URLs. Unlike the study/sample searches, this hits the read_run result domain, which is the only one that exposes downloadable read files (fastq_ftp, submitted_ftp, sra_ftp) plus checksums (fastq_md5), sizes (fastq_bytes), library_strategy and instrument_platform. Use ENA query syntax: 'tax_eq(TAXID)' for an exact species (e.g., tax_eq(1280) = Staphylococcus aureus), 'tax_tree(TAXID)' to include descendants, optionally combined with 'library_strategy="WGS"'. Returns each run with FTP URLs you can wget/curl directly. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENAPortal_search_runs_fastq',
     'arguments': {'query': '...'}},
)
```

## `ENAPortal_search_samples` {#ENAPortal_search_samples}

Searches ENA biological samples by taxonomy or text.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `query`
- **Keywords:** ena, biosample, samples, insdc
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | ENA search query. Examples: 'tax_tree(562)' (E. coli), 'description="liver"', 'tax_tree(9606) AND description="tumor"'. |
| `limit` | integer \| null | no | Maximum results to return (1-100, default 10). |
| `fields` | string \| null | no | Comma-separated fields. Default: 'sample_accession,sample_alias,description,tax_id,scientific_name,first_public'. |

??? note "Upstream description"

    Search the European Nucleotide Archive (ENA) for biological samples by taxonomy or text query. Returns sample accessions, descriptions, taxonomy, and publication dates. Useful for finding source material metadata across the world's nucleotide archives (INSDC). Example: 'tax_tree(562)' finds all E. coli samples; 'description="liver"' finds samples with liver in the description.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENAPortal_search_samples', 'arguments': {'query': '...'}},
)
```

## `ENAPortal_search_studies` {#ENAPortal_search_studies}

Searches ENA sequencing studies by text or taxonomy.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `query`
- **Keywords:** ena, studies, sequencing projects
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | ENA search query. IMPORTANT: Plain text keywords (e.g., 'cancer') will cause HTTP 400. Always wrap text keywords with description="...". Examples: 'description="cancer"', 'tax_tree(9606)' (human), 'description="transcriptomics" AND tax_tree(10090)' (mouse transcriptomics). |
| `limit` | integer \| null | no | Maximum results to return (1-100, default 10). |
| `fields` | string \| null | no | Comma-separated fields to return. Default: 'study_accession,study_title,center_name,first_public,description'. Other fields: 'last_updated,tax_id,scientific_name'. |

??? note "Upstream description"

    Search the European Nucleotide Archive (ENA) for sequencing studies using text queries or taxonomy filters. Returns study accessions, titles, submitting centers, and publication dates. ENA is the world's largest repository of nucleotide sequences (&gt;3 billion records). Supports powerful query syntax: use 'description="keyword"' for text search, 'tax_tree(TAXID)' for taxonomy-based queries, or combine with AND/OR. Example: searching 'description="SARS-CoV-2"' returns COVID-19 sequencing studies.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENAPortal_search_studies', 'arguments': {'query': '...'}},
)
```

## `ENCODE_get_biosample` {#ENCODE_get_biosample}

Gets ENCODE biosample metadata by accession (type, organism, donor, treatments).

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** encode, biosample, donor, cell line
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | ENCODE biosample accession identifier (format: ENCBS\#\#\#\#\#\#, e.g., 'ENCBS000AAA', 'ENCBS130ENC'). Find accessions using ENCODE_search_biosamples or from experiment metadata. |

??? note "Upstream description"

    Get detailed metadata for a specific ENCODE biosample by its accession ID. Returns comprehensive information including biosample type, organism, life stage, age, sex, tissue/cell line name, donor information, treatments, genetic modifications, culture conditions, and characterizations. Essential for understanding sample provenance, experimental conditions, and biological context of ENCODE datasets. Find biosample accessions using ENCODE_search_biosamples (format: ENCBS\#\#\#\#\#\#).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_get_biosample', 'arguments': {'accession': '...'}},
)
```

## `ENCODE_get_experiment` {#ENCODE_get_experiment}

Gets ENCODE experiment metadata by accession including assay, target, biosample, replicates and files.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** encode, experiment, chip-seq, assay, replicates
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | ENCODE experiment accession identifier (format: ENCSR\#\#\#\#\#\#, e.g., 'ENCSR000AKS', 'ENCSR000CAG'). Find accessions using ENCODE_search_experiments. |

??? note "Upstream description"

    Get comprehensive metadata for a specific ENCODE experiment by its accession ID. Returns detailed information including assay type, target protein, biosample details, replicates, files, quality metrics, analysis pipeline, lab, award, and publication references. Use this to understand experimental design, data quality, and available files before downloading. Find experiment accessions using ENCODE_search_experiments (format: ENCSR\#\#\#\#\#\#).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_get_experiment',
     'arguments': {'accession': '...'}},
)
```

## `ENCODE_get_file` {#ENCODE_get_file}

Gets ENCODE file metadata by accession including format, output type, assembly and download URL.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** encode, file, bam, bigwig, download
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | ENCODE file accession identifier (format: ENCFF\#\#\#\#\#\#, e.g., 'ENCFF001JXO'). Find accessions using ENCODE_list_files. |

??? note "Upstream description"

    Get detailed metadata for a specific ENCODE file by its accession ID. Returns comprehensive file information including file format, output type, dataset accession, biological/technical replicate numbers, file size, MD5 checksum, download URL (href), genome assembly, quality metrics, and processing status. Essential for accessing download URLs, verifying file integrity, and understanding file provenance. Find file accessions using ENCODE_list_files (format: ENCFF\#\#\#\#\#\#).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_get_file', 'arguments': {'accession': '...'}},
)
```

## `ENCODE_list_files` {#ENCODE_list_files}

Lists ENCODE data files filtered by format, output type or assay.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** encode, fastq, bam, bigwig, peaks
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `file_type` | string | no | File format filter (e.g., 'fastq', 'bam', 'bigWig', 'bed'). Common formats: fastq (raw reads), bam (alignments), bigWig (coverage), bed (peaks/regions). |
| `assay_title` | string | no | Assay filter (e.g., 'ChIP-seq', 'RNA-seq'). Filters files by the experimental assay that generated them. |
| `limit` | integer | no | Maximum number of results to return (1–100). Default `10` |

??? note "Upstream description"

    List ENCODE data files with filters by file format, output type, or assay. Returns file metadata including accessions, formats (FASTQ, BAM, bigWig), output types (alignments, peaks, signal), and file sizes. Use this to programmatically discover downloadable artifacts for specific experiments or data types. Find file accessions to use with ENCODE_get_file for detailed file information and download URLs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_list_files', 'arguments': {}},
)
```

## `ENCODE_search_biosamples` {#ENCODE_search_biosamples}

Searches ENCODE biosamples by organism, type or treatment.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** encode, biosample, cell line, tissue
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `organism` | string | no | Organism filter (e.g., 'Homo sapiens', 'Mus musculus'). |
| `biosample_type` | string | no | Biosample classification filter — maps to ENCODE's biosample_ontology.classification field (e.g., 'cell line', 'tissue', 'primary cell', 'in vitro differentiated cells'). |
| `treatment` | string | no | Treatment filter (e.g., 'interferon gamma', 'ethanol'). Use '*' to find any treated samples. |
| `limit` | integer | no | Maximum number of results to return (1–100). Default `10` |

??? note "Upstream description"

    Search ENCODE biosamples (cell lines, tissues, primary cells) by organism, biosample type, or treatment. Returns biosample accessions, types, organism, life stage, and descriptions. Use this to discover biological samples used in ENCODE experiments, find cell lines or tissues of interest, or identify treated vs. control samples. Find biosample accessions to use with ENCODE_get_biosample for detailed characterization data.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_search_biosamples', 'arguments': {}},
)
```

## `ENCODE_search_experiments` {#ENCODE_search_experiments}

Searches ENCODE experiments by assay, target, organism or status.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** encode, chip-seq, atac-seq, rna-seq, experiment search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `assay_title` | string | no | Assay name filter (e.g., 'TF ChIP-seq', 'Histone ChIP-seq', 'ATAC-seq', 'RNA-seq', 'Hi-C'). Leave empty to search all assays. Use 'TF ChIP-seq' for transcription factor binding experiments. |
| `target` | string | no | Target protein/factor filter (e.g., 'CTCF', 'H3K4me3', 'POLR2A'). Use for TF ChIP-seq experiments. |
| `organism` | string | no | Organism filter — maps to ENCODE's organism scientific name field (e.g., 'Homo sapiens', 'Mus musculus', 'Drosophila melanogaster'). |
| `status` | string | no | Record status filter. Use 'released' for public data (default), 'archived' for older data. Default `released` |
| `limit` | integer | no | Maximum number of results to return (1–100). Default `10` |

??? note "Upstream description"

    Search ENCODE functional genomics experiments (e.g., ChIP-seq, ATAC-seq, RNA-seq) by assay type, target protein, organism, or status. Returns experiment accessions, descriptions, and metadata. Use this to discover datasets for specific assays, transcription factors, or cell types. Find experiment accessions to use with ENCODE_get_experiment for detailed information or ENCODE_list_files for downloading data files.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_search_experiments', 'arguments': {}},
)
```

## `ENCODE_search_microrna_experiments` {#ENCODE_search_microrna_experiments}

Searches ENCODE microRNA-seq experiments.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** microrna-seq, mirna, encode, small rna
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `biosample_term_name` | string \| null | no | Biosample name filter (e.g., 'K562', 'GM12878', 'liver', 'heart'). Leave empty for all. |
| `organism` | string | no | Organism scientific name (e.g., 'Homo sapiens', 'Mus musculus'). Default `Homo sapiens` |
| `limit` | integer | no | Maximum number of results to return. Default `25` |

??? note "Upstream description"

    Search ENCODE microRNA-seq experiments. MicroRNA-seq profiles small regulatory RNAs (miRNAs, ~22nt) that post-transcriptionally regulate gene expression by targeting mRNA 3' UTRs. Returns experiment accessions, biosample details, and metadata. Use this to find miRNA expression profiling data for specific cell types or tissues. ENCODE has 400+ released microRNA-seq experiments across diverse cell types.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_search_microrna_experiments', 'arguments': {}},
)
```

## `ENCODE_search_rnaseq_experiments` {#ENCODE_search_rnaseq_experiments}

Searches ENCODE RNA-seq experiments by biosample and assay type.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** rna-seq, encode, polya, transcriptome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `biosample_term_name` | string \| null | no | Biosample name from ENCODE ontology (e.g., 'K562', 'HepG2', 'GM12878', 'liver', 'brain'). Leave empty to search all biosamples. |
| `biosample` | string \| null | no | Alias for biosample_term_name. Cell type or tissue (e.g., 'K562', 'liver'). |
| `cell_type` | string \| null | no | Alias for biosample_term_name. Cell type (e.g., 'K562', 'HepG2', 'GM12878'). |
| `tissue` | string \| null | no | Alias for biosample_term_name. Tissue type (e.g., 'liver', 'brain', 'heart'). |
| `organism` | string | no | Organism scientific name (e.g., 'Homo sapiens', 'Mus musculus'). Default `Homo sapiens` |
| `assay_type` | string | no | RNA-seq assay type: 'total RNA-seq' (default, all RNA), 'polyA plus RNA-seq' (mRNA), 'small RNA-seq', 'microRNA-seq'. Aliases: 'total', 'polya', 'mirna', 'small'. Default `total RNA-seq` |
| `limit` | integer | no | Maximum number of results to return (1-100). Default `25` |

??? note "Upstream description"

    Search ENCODE RNA-seq experiments by biosample (cell type or tissue), organism, and assay type (total RNA-seq, polyA plus RNA-seq, small RNA-seq, microRNA-seq). Returns experiment accessions, biosample summaries, lab names, and metadata from the ENCODE project. Use this to find ENCODE-curated RNA-seq data for specific cell lines or tissues.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ENCODE_search_rnaseq_experiments', 'arguments': {}},
)
```

## `FourDN_get_download_url` {#FourDN_get_download_url}

Gets the download URL and DRS endpoint for a 4DN data file by accession (requires 4DN access key).

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `file_accession`
- **Keywords:** 4dn, hi-c, drs, file download, 3d genome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | One of `download_file_url` Default `download_file_url` |
| `file_accession` | string | yes | 4DN file accession (e.g., '4DNFIIA7E3HL'). Obtain by searching with FourDN_search_data (item_type='File'). |

??? note "Upstream description"

    Get download URL and DRS (Data Repository Service) API endpoint for 4DN files. Prerequisites: (1) Requires file accession - obtain by using FourDN_search_data with item_type='File'. (2) Requires free 4DN account - create at data.4dnucleome.org, then generate access key in your profile settings. Returns download URL and instructions for command-line access with curl. Supports GA4GH DRS standard for interoperability. Use for: programmatic file download, batch processing, workflow integration.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FourDN_get_download_url',
     'arguments': {'file_accession': '...'}},
)
```

## `FourDN_get_experiment_metadata` {#FourDN_get_experiment_metadata}

Gets 4DN experiment metadata: assay type (Hi-C, ATAC-seq, ChIP-seq), biosource, protocol and output files.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `experiment_accession`
- **Keywords:** 4dn, hi-c, experiment metadata, chip-seq, atac-seq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | One of `get_experiment_metadata` Default `get_experiment_metadata` |
| `experiment_accession` | string | yes | 4DN experiment accession (e.g., '4DNEXHVF8WA9'). Obtain by searching with FourDN_search_data (item_type='Experiment'). |
| `include_full_metadata` | boolean | no | Include complete API response in 'metadata' field (default: false). Set to true only if you need all raw fields. Omitting this reduces response size significantly. Default `False` |

??? note "Upstream description"

    Get metadata for 4DN experiments including experimental design, biosource information, protocol details, and associated files. Returns experiment type (Hi-C variants, ATAC-seq, ChIP-seq), cell type, processing status, and list of output files. Prerequisites: Requires experiment accession - obtain by using FourDN_search_data with item_type='Experiment'. Use for: understanding experimental context, finding all files from experiment, checking data provenance.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'FourDN_get_experiment_metadata',
     'arguments': {'experiment_accession': '...'}},
)
```

## `GDC_list_files` {#GDC_list_files}

Lists GDC files filtered by data type (e.g. expression quantification) for analysis pipelines.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** gdc, tcga files, data type, download, genomic data
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `data_type` | string | no | Data type filter (e.g., 'Gene Expression Quantification'). |
| `size` | integer | no | Number of results (1–100). Default `10` |
| `offset` | integer | no | Offset for pagination (0-based). Default `0` |

??? note "Upstream description"

    List GDC files filtered by data_type and other fields. Use to identify downloadable artifacts (e.g., expression quantification) for analysis pipelines.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GDC_list_files', 'arguments': {}},
)
```

## `GDC_list_projects` {#GDC_list_projects}

Lists GDC projects (TCGA, TARGET, CGCI) with case and file counts, filterable by program.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** tcga, target, cancer genomics, cohorts, gdc
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `program` | string | no | Filter by program name (e.g., 'TCGA', 'TARGET') |
| `size` | integer | no | Number of results (1–100) Default `20` |

??? note "Upstream description"

    List GDC projects (TCGA, TARGET, etc.) with case/file counts. Use to discover available cancer cohorts and their data volume. Filter by program (TCGA, TARGET, CGCI, etc.).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GDC_list_projects', 'arguments': {}},
)
```

## `GEO_get_dataset_details` {#GEO_get_dataset_details}

Gets metadata for a GEO series (GSE) including platform, samples and supplementary links.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `geo_id`
- **Keywords:** geo, gse, series, platform
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `geo_id` | string | yes | GEO dataset accession (numeric part only, e.g., '200291249' for GSE291249). Get IDs from GEO search results. |

??? note "Upstream description"

    Get detailed metadata for a specific GEO dataset by its accession (GSE ID). Returns comprehensive information including title, summary, experiment type, platform, organism, sample count, publication date, and supplementary data links. Use this after searching with GEO_search_methylation_datasets or GEO_search_chipseq_datasets to get full details about a dataset of interest.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GEO_get_dataset_details', 'arguments': {'geo_id': '...'}},
)
```

## `geo_get_dataset_info` {#geo_get_dataset_info}

Gets title, summary and metadata for a GEO dataset.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `dataset_id`
- **Keywords:** geo, gse, gds, dataset metadata, functional genomics
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `dataset_id` | string | yes | GEO dataset ID (e.g., 'GDS1234', 'GSE12345') |

??? note "Upstream description"

    Get detailed information about a specific GEO dataset including title, summary, and metadata.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'geo_get_dataset_info',
     'arguments': {'dataset_id': '...'}},
)
```

## `geo_get_sample_info` {#geo_get_sample_info}

Gets sample characteristics and experimental conditions for a GEO dataset.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `dataset_id`
- **Keywords:** geo, gsm, sample metadata, experimental conditions
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `dataset_id` | string | yes | GEO dataset ID (e.g., 'GDS1234', 'GSE12345') |

??? note "Upstream description"

    Get sample information for a GEO dataset including sample characteristics and experimental conditions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'geo_get_sample_info', 'arguments': {'dataset_id': '...'}},
)
```

## `geo_list_supplementary_files` {#geo_list_supplementary_files}

Lists downloadable supplementary/raw files of a GEO GSE or GSM with size, type and download URL.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** geo, supplementary files, raw data, ftp, download
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | GEO Series (GSE...) or Sample (GSM...) accession. Examples: 'GSE42657', 'GSE1000', 'GSM1045442'. |

??? note "Upstream description"

    List the actual downloadable supplementary/raw files of a GEO Series (GSE) or Sample (GSM), with per-file name, byte size, type, modification time, and a direct download URL. GEO's E-utilities esummary only reports a type code and the base FTP directory, never the per-file manifest. For a Series this reads the structured suppl/filelist.txt TSV (Archive/File, Name, Time, Size, Type); for a Sample (no filelist.txt) it parses the suppl/ directory listing. Example: GSE42657 returns GSE42657_RAW.tar (TAR, 71680 bytes) and GPL8179_humanMI_V2_R0_XS0000124-MAP.txt.gz (TXT). Keyless; no NCBI API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'geo_list_supplementary_files',
     'arguments': {'accession': '...'}},
)
```

## `geo_search_datasets` {#geo_search_datasets}

Searches GEO datasets by query, organism, study type or platform and returns IDs and basic info.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** geo, dataset search, microarray, rna-seq, platform
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | no | Search query terms (e.g., 'cancer', 'diabetes', 'microarray') |
| `organism` | string | no | Organism name (e.g., 'Homo sapiens', 'Mus musculus') Default `` |
| `study_type` | string | no | Type of study (e.g., 'Expression profiling by array', 'Expression profiling by high throughput sequencing') Default `` |
| `platform` | string | no | Platform used (e.g., 'GPL570', 'GPL96') Default `` |
| `limit` | integer | no | Maximum number of results to return Default `50` |

??? note "Upstream description"

    Search GEO datasets by query terms, organism, study type, or platform. Returns dataset IDs and basic information.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'geo_search_datasets', 'arguments': {}},
)
```

## `GEO_search_rnaseq_datasets` {#GEO_search_rnaseq_datasets}

Searches NCBI GEO for RNA-seq datasets by keyword, tissue or gene.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** geo, rna-seq, transcriptomics, gse
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string \| null | no | Search query (e.g., 'BRCA1 breast cancer', 'T cell activation', 'Alzheimer disease'). Leave empty to search broadly. |
| `organism` | string | no | Organism scientific name (e.g., 'Homo sapiens', 'Mus musculus'). Default `Homo sapiens` |
| `limit` | integer | no | Maximum number of results to return (1-100). Default `20` |

??? note "Upstream description"

    Search NCBI GEO for RNA-seq gene expression datasets by query keyword, disease, tissue, or gene of interest. Returns dataset accessions (GSE IDs), titles, summaries, organism, and sample counts. Use this to find transcriptomics experiments relevant to a biological question. Covers bulk RNA-seq and total RNA-seq studies across all tissues and diseases.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GEO_search_rnaseq_datasets', 'arguments': {}},
)
```

## `GSA_get_accession` {#GSA_get_accession}

Looks up a China NGDC Genome Sequence Archive CRA accession returning title and linked BioProject/BioSample.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** gsa, ngdc, cra, sequencing archive, bioproject
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | GSA accession, e.g. 'CRA002926'. |

??? note "Upstream description"

    Look up a GSA (Genome Sequence Archive) accession by its CRA-format ID. GSA, run by China's National Genomics Data Center, is a major raw-sequencing archive alongside NCBI SRA/EBI ENA/DDBJ for a large and growing share of Chinese-origin genomic datasets. Returns title, linked BioProject accession, associated publication (title/journal/DOI/PubMed ID), file count/size, and direct HTTPS/FTP download URLs. Example: accession='CRA002926' returns a gynecological-mucosal-melanoma sequencing dataset linked to its Frontiers in Oncology publication.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GSA_get_accession', 'arguments': {'accession': '...'}},
)
```

## `GTEx_get_sample_info` {#GTEx_get_sample_info}

Gets GTEx sample and subject metadata: tissue, donor age/sex, ischemic time, RIN and QC metrics.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** gtex, sample metadata, rin, donor, qc
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `get_sample_info` Default `get_sample_info` |
| `sample_id` | array of string | no | Optional: GTEx sample IDs |
| `subject_id` | array of string | no | Optional: GTEx subject/donor IDs |
| `tissue_site_detail_id` | array of string | no | Optional: Filter by tissues |
| `sex` | string | no | Optional: Filter by donor sex One of `male`, `female` |
| `age_bracket` | array of string | no | Optional: Filter by age brackets |
| `dataset_id` | string | no | One of `gtex_v8`, `gtex_v10`, `gtex_snrnaseq_pilot` Default `gtex_v8` |
| `page` | integer | no | Default `0` |
| `items_per_page` | integer | no | Default `250` |

??? note "Upstream description"

    Get detailed GTEx sample and subject metadata. Returns sample IDs, tissue types, donor demographics (age, sex), ischemic time, RIN scores, Hardy scale, pathology notes, and QC metrics. Filter by sample ID, subject ID, tissue, demographics, or data type. Use for: sample selection, cohort design, QC filtering, demographic analysis, understanding data provenance. Essential for selecting appropriate samples for analysis.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'GTEx_get_sample_info', 'arguments': {}},
)
```

## `hca_get_file_manifest` {#hca_get_file_manifest}

Retrieves downloadable file list for a Human Cell Atlas project by project ID.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `project_id`
- **Keywords:** human cell atlas, single-cell, file manifest, download
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `action` | string | no | The specific action to perform. Must be set to 'get_file_manifest'. One of `get_file_manifest` Default `get_file_manifest` |
| `project_id` | string | yes | The unique UUID of the project (entryId) for which to retrieve files. This ID is typically obtained from the 'hca_search_projects' tool. |
| `limit` | integer | no | The maximum number of files to list. The default is 10. Default `10` |

??? note "Upstream description"

    Retrieves a list of downloadable files for a specific Human Cell Atlas (HCA) project identified by its Project ID.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'hca_get_file_manifest',
     'arguments': {'project_id': '...'}},
)
```

## `hca_search_projects` {#hca_search_projects}

Searches Human Cell Atlas DCP single-cell projects by organ and disease.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** human cell atlas, single-cell, projects, organ
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `action` | string | no | The specific action to perform. Must be set to 'search_projects'. One of `search_projects` Default `search_projects` |
| `organ` | string | no | The organ to filter projects by. Examples include 'heart', 'liver', 'brain', 'kidney'. This field is optional. |
| `disease` | string | no | The disease state to filter by. Examples include 'normal', 'cancer', 'covid-19'. This field is optional. |
| `limit` | integer | no | The maximum number of projects to return. The default is 10. Use this to control the size of the result set. Default `10` |

??? note "Upstream description"

    Searches for single-cell projects in the Human Cell Atlas (HCA) Data Coordination Platform (DCP). This tool allows filtering projects based on specific organs and disease states.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'hca_search_projects', 'arguments': {}},
)
```

## `IGSR_list_data_collections` {#IGSR_list_data_collections}

Lists IGSR data collections (1000 Genomes phases, 30x resequencing, HGSVC, HPRC).

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** 1000 genomes, igsr, data collections, pangenome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type (fixed: list_data_collections). One of `list_data_collections` Default `list_data_collections` |
| `limit` | integer | no | Maximum number of results to return. Default `50` |

??? note "Upstream description"

    List available data collections from the International Genome Sample Resource (IGSR). Data collections include the original 1000 Genomes Project phases, the 30x high-coverage resequencing, Human Genome Structural Variation Consortium (HGSVC), Human Pangenome Reference Consortium (HPRC), and other related studies. Returns collection codes, titles, sample/population counts, and available data types (sequence, alignment, variants). Use this to discover what genomics datasets are available through IGSR.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IGSR_list_data_collections', 'arguments': {}},
)
```

## `IGSR_search_samples` {#IGSR_search_samples}

Searches 1000 Genomes samples by population, data collection or name with sex and BioSample IDs.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** 1000 genomes, samples, na12878, biosample
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type (fixed: search_samples). One of `search_samples` Default `search_samples` |
| `population` | string \| null | no | Filter by population code (e.g., 'CHB', 'YRI', 'GBR', 'CEU'). Use IGSR_search_populations to find valid codes. |
| `data_collection` | string \| null | no | Filter by data collection title (e.g., '1000 Genomes 30x on GRCh38', '1000 Genomes on GRCh38'). |
| `sample_name` | string \| null | no | Look up specific sample by name (e.g., 'HG00096', 'NA12878'). |
| `limit` | integer | no | Maximum number of results to return. Default `25` |

??? note "Upstream description"

    Search samples from the 1000 Genomes Project by population, data collection, or sample name. Returns sample names (e.g., HG00096, NA12878), sex, BioSample IDs, population assignments, and available data collections. Use this to find samples from specific populations for variant analysis, to identify samples with particular sequencing data types (WGS, exome, RNA-seq), or to look up specific well-characterized samples like NA12878 (reference genome sample).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IGSR_search_samples', 'arguments': {}},
)
```

## `IntOGen_list_cohorts` {#IntOGen_list_cohorts}

Lists IntOGen cancer cohorts with data source, sample counts, tumor type and age group, optionally filtered by cancer type.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** cancer cohorts, tcga, icgc, hartwig, sample counts
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `cancer_type` | string \| null | no | Optional cancer type code to filter cohorts. If omitted, returns all 271 cohorts. Examples: 'BRCA', 'LUAD', 'COADREAD'. |

??? note "Upstream description"

    List available cancer cohorts from IntOGen with sample counts and tumor types. Returns cohort IDs, names, data sources (TCGA, ICGC, Hartwig, etc.), sample counts, tumor type (Primary/Metastasis), and age group (Adult/Pediatric). Optionally filter by cancer type. IntOGen includes 271 cohorts from major genomic studies. Example: cancer_type='BRCA' returns 14 breast cancer cohorts from TCGA, Hartwig, ICGC, and others.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'IntOGen_list_cohorts', 'arguments': {}},
)
```

## `MGnify_get_samples` {#MGnify_get_samples}

Lists MGnify samples with BioSample ID, location, host and environment metadata.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** biosample, metagenomic samples, envo, geolocation
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sample_accession` | string \| null | no | Fetch a single sample by its MGnify/ENA accession (e.g., 'SRS10016989'). If provided, study_accession/biome/page are ignored. |
| `study_accession` | string \| null | no | Filter samples belonging to a study (e.g., 'MGYS00002008'). |
| `biome` | string \| null | no | Filter by biome name lineage (e.g., 'root:Host-associated:Human:Digestive system'). |
| `page` | integer \| null | no | Page number (default 1). |
| `page_size` | integer \| null | no | Samples per page (default 25, max 100). |

??? note "Upstream description"

    List biological SAMPLES in MGnify (EBI Metagenomics) with full geographic/host/environment provenance. Returns each sample's accession, BioSample id, latitude, longitude, geo-loc-name, collection-date, host-tax-id, species, and ENVO environment biome/feature/material - the sample-level metadata table that the study, analysis, genome, and biome tools never expose. Pass a single 'sample_accession' (e.g. 'SRS10016989') to fetch one sample, or omit it to page through samples optionally filtered by 'study_accession' or 'biome'. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MGnify_get_samples', 'arguments': {}},
)
```

## `MGnify_get_study_detail` {#MGnify_get_study_detail}

Returns MGnify study details with abstract, BioProject, sample count and biomes.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `study_accession`
- **Keywords:** mgnify study, bioproject, metagenomics
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `study_accession` | string | yes | MGnify study accession. Examples: 'MGYS00002008', 'MGYS00005292'. |

??? note "Upstream description"

    Get detailed information about a specific MGnify metagenomics study. Returns study name, abstract, associated BioProject, centre, public/private status, sample count, and biomes. Complementary to MGnify_search_studies which lists studies; use MGnify_list_analyses / MGnify_list_analysis_downloads for exact analyses/downloads counts (this endpoint does not provide them). Example: MGYS00002008 is a Tara Oceans marine metagenomics study.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MGnify_get_study_detail',
     'arguments': {'study_accession': '...'}},
)
```

## `MGnify_get_taxonomy` {#MGnify_get_taxonomy}

Returns taxonomic composition and abundances from a MGnify analysis using SSU or LSU rRNA markers.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `analysis_id`
- **Keywords:** 16s rrna, taxonomic profile, abundance, lsu
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `analysis_id` | string | yes | MGnify analysis accession (e.g., 'MGYA00585482'). Find IDs via MGnify_list_analyses. |
| `rna_type` | string | no | rRNA marker type: 'ssu' for 16S/18S (default) or 'lsu' for 23S/28S. Default `ssu` |
| `page_size` | integer | no | Number of taxa to return per page (default 25, max 100). Default `25` |

??? note "Upstream description"

    Get taxonomic composition from a MGnify metagenomics analysis. Returns organisms detected with their abundance counts, taxonomic lineage (domain/phylum/class/order/family/genus/species), and rank. Supports SSU (16S/18S) and LSU (23S/28S) rRNA marker genes. Use MGnify_list_analyses to find analysis IDs for a study. Example: MGYA00585482 returns Bacteria, Proteobacteria, Firmicutes lineages with counts.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MGnify_get_taxonomy',
     'arguments': {'analysis_id': '...'}},
)
```

## `MGnify_list_analyses` {#MGnify_list_analyses}

Lists MGnify analyses for a study accession with taxonomic and functional outputs.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `study_accession`
- **Keywords:** metagenomics, mgnify analyses, study accession
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `study_accession` | string | yes | MGnify study accession (e.g., 'MGYS00000001'). |
| `size` | integer | no | Number of records per page (1–100). Default `10` |

??? note "Upstream description"

    List analyses associated with a study accession (taxonomic/functional outputs). Use to enumerate available processed results for programmatic retrieval.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MGnify_list_analyses',
     'arguments': {'study_accession': '...'}},
)
```

## `MGnify_list_analysis_downloads` {#MGnify_list_analysis_downloads}

Lists MGnify analysis output files such as predicted CDS and taxonomy tables with download URLs.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `analysis_id`
- **Keywords:** download files, predicted cds, fasta, mgnify
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `analysis_id` | string | yes | MGnify analysis accession (e.g., 'MGYA00585482'). Find IDs via MGnify_list_analyses. |
| `page_size` | integer \| null | no | Number of download entries to return (default 100, max 100). |

??? note "Upstream description"

    List the downloadable result files produced by a MGnify analysis. Returns each output file's id, label/description, file format, compression, group type, and a direct download URL - covering taxonomy and functional outputs such as predicted CDS (.faa), predicted ORF (.ffn), processed contigs (.fasta), combined eggNOG/InterPro/antiSMASH annotation GFF, eggNOG (emapper) TSV, and Diamond TSV. The parsed-annotation tools (taxonomy, GO, InterPro) summarize content but cannot enumerate or point to the actual files; this tool does. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MGnify_list_analysis_downloads',
     'arguments': {'analysis_id': '...'}},
)
```

## `MGnify_search_studies` {#MGnify_search_studies}

Searches MGnify metagenomics/microbiome studies by biome or keyword.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** metagenomics, microbiome, biome, bioproject
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `biome` | string | no | Biome identifier (e.g., 'root:Host-associated'). |
| `search` | string | no | Keyword to search in study titles/descriptions. |
| `size` | integer | no | Number of records per page (1–100). Default `10` |

??? note "Upstream description"

    Search MGnify metagenomics/microbiome studies by biome or keyword. Returns study accessions (MGYS...), names, bioproject IDs, biome names, and analysis counts. Use to find datasets from specific environments (e.g., gut microbiome, soil, ocean) before retrieving analyses or taxonomy. Example: biome='root:Host-associated' returns gut, skin, and oral studies.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MGnify_search_studies', 'arguments': {}},
)
```

## `MPD_get_phenotype_data` {#MPD_get_phenotype_data}

Searches ENCODE experiment records that mention a given mouse strain.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** encode experiments, mouse strain, phenome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `strain` | string | no | Mouse strain (e.g., C57BL/6J, BALB/c, DBA/2J) Default `C57BL/6J` |
| `phenotype_category` | string | no | Not currently applied to the query -- ENCODE (this tool's data source) has no phenotype-category concept. Retained for forward compatibility with a future genuine MPD-backed implementation. Default `behavior` |
| `limit` | integer | no | Number of results to return Default `10` |

??? note "Upstream description"

    Search ENCODE experiment records mentioning a given mouse strain (used as a Mouse Phenome Database alternative -- MPD's own REST API has no simple strain+category phenotype search; see https://phenome.jax.org/api for genuine curated MPD phenotype data). 'phenotype_category' is not applied -- ENCODE has no equivalent concept.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'MPD_get_phenotype_data', 'arguments': {}},
)
```

## `NCBI_SRA_get_download_urls` {#NCBI_SRA_get_download_urls}

Returns FTP, AWS S3 and web URLs for downloading SRA run data.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accessions`
- **Keywords:** sra, fastq, download, aws s3
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | any | no | Operation type (fixed: get_download_urls) |
| `accessions` | any | yes | SRA run accession(s) (e.g., 'SRR000001', 'ERR000001', 'DRR000001'). Can be single accession or array of accessions. |

??? note "Upstream description"

    Get FTP and cloud storage URLs for downloading FASTQ sequencing data from SRA runs. Returns multiple access methods: NCBI FTP URLs for direct .sra file download, AWS S3 URLs for cloud-based access (faster for large-scale analysis), and NCBI web URLs for browser access. Note: Downloaded .sra files require SRA Toolkit (fastq-dump or fasterq-dump) to convert to standard FASTQ format. Essential for programmatic access to petabytes of public sequencing data. Works with SRR (NCBI), ERR (ENA), and DRR (DDBJ) accessions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBI_SRA_get_download_urls',
     'arguments': {'accessions': '...'}},
)
```

## `NCBI_SRA_get_run_info` {#NCBI_SRA_get_run_info}

Returns SRA run metadata including platform, instrument, library strategy, layout and sample details.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accessions`
- **Keywords:** sra run, library strategy, illumina, wgs, rna-seq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | any | no | Operation type (fixed: get_run_info) |
| `accessions` | any | yes | SRA run accession(s) (e.g., 'SRR000001', 'ERR000001', 'DRR000001'). Can be single accession or array of accessions. Use SRR for NCBI, ERR for ENA, DRR for DDBJ. |

??? note "Upstream description"

    Get detailed metadata for SRA run accessions including sequencing platform (ILLUMINA, Oxford Nanopore), instrument model, library strategy (RNA-Seq, WGS), library source (genomic, transcriptomic), library layout (paired/single), sample information, study details, organism, and data statistics (total spots, bases). Essential for understanding experimental design and data characteristics before downloading large FASTQ files. Use SRA run accessions (SRR000001, ERR000001, DRR000001) from NCBI_SRA_search_runs or known accessions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBI_SRA_get_run_info',
     'arguments': {'accessions': '...'}},
)
```

## `NCBI_SRA_link_to_biosample` {#NCBI_SRA_link_to_biosample}

Links SRA runs to BioSample records with sample attributes, tissue, cell line and donor information.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accessions`
- **Keywords:** biosample, sra, sample metadata
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | any | no | Operation type (fixed: link_to_biosample) |
| `accessions` | any | yes | SRA UID(s) from NCBI_SRA_search_runs. Can be single UID or array of UIDs. These are NCBI internal IDs (numeric), not SRR accessions. |

??? note "Upstream description"

    Link SRA sequencing runs to their associated BioSample records to access biological sample metadata including sample attributes, tissue type, cell line, treatment conditions, donor information, and collection details. BioSample provides rich contextual information about the biological specimens used in sequencing experiments. Essential for understanding experimental context and sample provenance. Returns BioSample UIDs which can be used with NCBI BioSample API to retrieve detailed sample metadata. Works with SRA UIDs from NCBI_SRA_search_runs.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBI_SRA_link_to_biosample',
     'arguments': {'accessions': '...'}},
)
```

## `NCBI_SRA_locate_run_files` {#NCBI_SRA_locate_run_files}

Resolves cloud download locations with file size and md5 for an SRA run via the SRA Data Locator.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accessions`
- **Keywords:** sra data locator, md5, cloud download, run files
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | any | no | Operation type (fixed: locate_run_files) |
| `accessions` | any | yes | SRA run accession(s) (e.g., 'SRR390728', 'ERR000001', 'DRR000001'). Single accession or array. Use SRR for NCBI, ERR for ENA, DRR for DDBJ. |

??? note "Upstream description"

    Resolve the current, verified cloud download location(s) plus authoritative file size and md5 checksum for an SRA run via the NCBI SRA Data Locator (SDL) service. Unlike NCBI_SRA_get_download_urls (which builds a guessed legacy ftp-trace ByRun path that now returns HTTP 404, and an S3 path with no size/md5/region), this returns the real object: https S3/GCP link, exact byte size, md5, storage region, and modification date. Use SRR/ERR/DRR run accessions. Example: SRR390728 returns size=195174182, md5=29a6a1a0dd0702f45225f2eb93c958b5, location s3:us-east-1 https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR390728/SRR390728. Keyless; no NCBI API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBI_SRA_locate_run_files',
     'arguments': {'accessions': '...'}},
)
```

## `NCBI_SRA_search_runs` {#NCBI_SRA_search_runs}

Searches NCBI SRA runs by study, organism, strategy, platform or library source.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** sra search, ngs, chip-seq, oxford nanopore
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | any | no | Operation type (fixed: search) |
| `study` | string | no | Study accession to search (e.g., 'SRP000001', 'ERP000001', 'DRP000001'). Searches in \[Study\] field. |
| `organism` | string | no | Organism name to filter (e.g., 'Homo sapiens', 'Mus musculus', 'SARS-CoV-2'). Searches in \[Organism\] field. |
| `strategy` | string | no | Sequencing strategy to filter. Common values: 'RNA-Seq' (transcriptomics), 'WGS' (whole genome), 'WXS' (exome), 'ChIP-Seq' (chromatin), 'ATAC-Seq' (accessibility), 'Bisulfite-Seq' (methylation), 'AMPLICON' (targeted), 'Hi-C' (3D genome) |
| `platform` | string | no | Sequencing platform to filter. Common values: 'ILLUMINA', 'OXFORD_NANOPORE', 'PACBIO_SMRT', 'ION_TORRENT', 'BGISEQ' |
| `source` | string | no | Library source to filter. Common values: 'GENOMIC' (DNA), 'TRANSCRIPTOMIC' (RNA), 'METAGENOMIC', 'METATRANSCRIPTOMIC', 'SYNTHETIC', 'VIRAL_RNA' |
| `query` | string | no | Free-form search query if not using structured filters. Use NCBI query syntax (e.g., '(Homo sapiens\[Organism\]) AND (RNA-Seq\[Strategy\])') |
| `limit` | integer | no | Maximum number of results to return (default: 20, max: 100) Default `20` |
| `sort` | string | no | Sort order for results (default: relevance) One of `relevance`, `pub_date` Default `relevance` |

??? note "Upstream description"

    Search NCBI Sequence Read Archive (SRA) for NGS/RNA-seq sequencing runs by study accession, organism, sequencing strategy (RNA-Seq, WGS, ChIP-Seq), platform (ILLUMINA, OXFORD_NANOPORE), or library source. Returns SRA UIDs which can be used with other SRA tools to retrieve metadata and download links. Essential for discovering high-throughput sequencing data from public repositories. Use this to find relevant sequencing experiments before accessing the actual FASTQ data files.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NCBI_SRA_search_runs', 'arguments': {}},
)
```

## `Nextstrain_get_dataset` {#Nextstrain_get_dataset}

Returns metadata and phylogenetic summary for a Nextstrain pathogen dataset, including sequence counts.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `dataset`
- **Keywords:** nextstrain, phylogenetics, pathogen genomes, genomic epidemiology
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `dataset` | string | yes | Nextstrain dataset path. Examples: 'zika', 'ebola', 'flu/seasonal/h3n2/ha/2y', 'dengue/denv1', 'mpox/clade-IIb'. |

??? note "Upstream description"

    Get metadata and phylogenetic summary for a Nextstrain pathogen dataset. Returns title, last update date, number of sequences, data provenance, maintainers, root node attributes, and available coloring options. Useful for understanding the scope and currency of pathogen genomic surveillance data. Example: 'zika' returns 1044 sequences tracking Zika virus evolution, last updated Dec 2025.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Nextstrain_get_dataset', 'arguments': {'dataset': '...'}},
)
```

## `Nextstrain_list_datasets` {#Nextstrain_list_datasets}

Lists Nextstrain pathogen phylogenetic datasets grouped by pathogen.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** nextstrain, sars-cov-2, influenza, phylogeny
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `pathogen` | string \| null | no | Optional pathogen name filter. Examples: 'flu', 'ebola', 'zika', 'dengue', 'mpox'. Leave empty for all pathogens. |
| `datasets_per_pathogen` | integer \| null | no | Maximum number of dataset paths listed per pathogen (default 10). Set to 0 to list every dataset with no cap, or to a larger number (e.g. 100) to cover the biggest pathogens such as ncov and avian-flu. |

??? note "Upstream description"

    List available pathogen phylogenetic datasets from Nextstrain. Returns datasets grouped by pathogen (e.g., avian-flu, dengue, ebola, flu, measles, mpox, SARS-CoV-2, zika), each with its true 'dataset_count'. By default only the first 10 dataset paths per pathogen are listed; when that hides entries the response sets 'truncated': true with a 'truncation_note', and 'metadata.total_datasets' always reports the full catalogue size (not the number listed). Pass datasets_per_pathogen=0 to list every dataset. Optionally filter by pathogen name. Nextstrain provides real-time molecular epidemiology for public health, tracking how pathogens evolve and spread worldwide.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Nextstrain_list_datasets', 'arguments': {}},
)
```

## `OSDR_get_study` {#OSDR_get_study}

Returns metadata for a NASA OSDR space biology study including organism, mission, factors and assay type.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `study_id`
- **Keywords:** osdr, genelab, spaceflight, transcriptomics
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `study_id` | string | yes | OSDR study accession, e.g. 'OSD-1' or '1'. |

??? note "Upstream description"

    Retrieve full metadata for one NASA OSDR space biology study: title, description, organism, spaceflight mission, experimental factors, assay type, and publication details. Example: study_id='OSD-1' (Drosophila immune response aboard the Space Shuttle). Use OSDR_search_studies to find a study_id.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OSDR_get_study', 'arguments': {'study_id': '...'}},
)
```

## `OSDR_list_files` {#OSDR_list_files}

Lists downloadable raw and processed data files for a NASA OSDR study.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `study_id`
- **Keywords:** osdr files, genelab, microarray, download
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `study_id` | string | yes | OSDR study accession, e.g. 'OSD-1' or '1'. |
| `limit` | integer \| null | no | Max files to return, 1-200. Default 50. |

??? note "Upstream description"

    List the downloadable data files for one NASA OSDR study: raw and processed files across assay types, with category, size, and download path. Example: study_id='OSD-1' returns 73 files including microarray differential expression tables. Use OSDR_search_studies to find a study_id.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OSDR_list_files', 'arguments': {'study_id': '...'}},
)
```

## `OSDR_search_studies` {#OSDR_search_studies}

Searches NASA OSDR space biology studies across transcriptomics, proteomics and other assays.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `query`
- **Keywords:** osdr search, genelab, microgravity, omics
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Free-text search, e.g. 'microgravity bone loss', 'radiation Drosophila'. |
| `organism` | string \| null | no | Optional organism filter, e.g. 'Mus musculus', 'Homo sapiens', 'Arabidopsis thaliana'. |
| `limit` | integer \| null | no | Max studies to return, 1-100. Default 25. |

??? note "Upstream description"

    Full-text search NASA's Open Science Data Repository (OSDR, formerly GeneLab): ~1,000 space biology studies on how spaceflight and microgravity affect organisms, across transcriptomics, proteomics, microscopy, and physiology assays. Optionally filter by organism. Example: query='microgravity bone loss', organism='Mus musculus'. Results also include cross-referenced external datasets (e.g. GEO); only rows with usable_with_get_study true work with OSDR_get_study and OSDR_list_files.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'OSDR_search_studies', 'arguments': {'query': '...'}},
)
```

## `Pathoplexus_get_sequence_details` {#Pathoplexus_get_sequence_details}

Lists per-sequence metadata records (accession, country, date, host, lineage) for a Pathoplexus organism.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `organism`
- **Keywords:** accession, metadata, pathogen genome
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `organism` | string | yes | Pathoplexus organism slug: 'west-nile', 'ebola-zaire', 'ebola-sudan', 'cchf', 'mpox'. |
| `country` | string \| null | no | Filter by country (geoLocCountry), e.g. 'USA'. |
| `lineage` | string \| null | no | Filter by lineage. |
| `fields` | string \| null | no | Comma-separated metadata fields to return, e.g. 'accession,geoLocCountry,sampleCollectionDate,lineage'. Omit to return all available fields. |
| `limit` | integer \| null | no | Max sequence rows (default 50, max 1000). |
| `offset` | integer \| null | no | Row offset for pagination (default 0). |

??? note "Upstream description"

    List the individual sequence records (not aggregated counts) behind a Pathoplexus organism, from the LAPIS /details endpoint. Returns a per-sequence metadata table — accession, country (geoLocCountry), collection date (sampleCollectionDate), host, lineage, clade and any other requested fields — so you can enumerate the actual surveilled genomes rather than only their totals. Optionally filter by country and lineage and select which fields to return. Known organisms: 'west-nile', 'ebola-zaire', 'ebola-sudan', 'cchf', 'mpox'. No API key required.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Pathoplexus_get_sequence_details',
     'arguments': {'organism': '...'}},
)
```

## `PharmacoDB_get_molecular_profiling` {#PharmacoDB_get_molecular_profiling}

Lists which molecular profiling data types (RNA-seq, mutation, CNV) exist for a cell line and in which dataset.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** ccle, rna-seq, cnv, molecular profiling
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `operation` | string | no | Operation type One of `get_molecular_profiling` Default `get_molecular_profiling` |
| `cell_line_name` | string \| null | no | Cancer cell line name (e.g., 'MCF-7', 'A549'). Mutually exclusive with cell_line_id. |
| `cell_line_id` | integer \| null | no | PharmacoDB cell line database ID (e.g., 273 for MCF-7). Mutually exclusive with cell_line_name. |

??? note "Upstream description"

    Get the per-cell-line molecular profiling inventory from PharmacoDB: which molecular data types (e.g., Kallisto RNA-seq isoforms, rnaseq, rna microarray, mutation, cnv) are available for a given cancer cell line, and in which pharmacogenomics dataset (CCLE, CTRPv2, FIMM, GDSC1/2, gCSI, etc.), with the number of profiles available (num_prof). Use this to confirm a cell line has the omics layer you need before pulling biomarker associations. Look up by cell_line_name (e.g., 'MCF-7') or cell_line_id. Example: cell_line_name='MCF-7' (cell id 273) returns records like \{dataset CCLE, mDataType 'Kallisto_0.46.1.rnaseq', num_prof 1\} plus mutation/cnv/rna layers.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'PharmacoDB_get_molecular_profiling', 'arguments': {}},
)
```

## `Progenetix_list_cohorts` {#Progenetix_list_cohorts}

Lists Progenetix cancer CNV cohorts (e.g. TCGA) with IDs, sizes and data types.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** cnv, cohort, tcga, cancer genomics
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `limit` | integer | no | Maximum number of cohorts to return (default: 10). Default `10` |

??? note "Upstream description"

    List available cohorts in the Progenetix cancer CNV database. Returns cohort IDs, names, sizes, and data types. Progenetix organizes samples into cohorts like TCGA, DIPG studies, and curated cancer sets. Example: returns 6 cohorts including TCGA (517 samples) and DIPG (16 samples).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Progenetix_list_cohorts', 'arguments': {}},
)
```

## `Progenetix_search_biosamples` {#Progenetix_search_biosamples}

Searches Progenetix tumor biosamples by NCIt disease code, returning IDs, histology and stage.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `filters`
- **Keywords:** biosample, cnv profile, ncit, tumor
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `filters` | string | yes | NCIt ontology code to filter biosamples by cancer type. Examples: 'NCIT:C3058' (Glioblastoma), 'NCIT:C4017' (Breast Carcinoma), 'NCIT:C9145' (Acute Myeloid Leukemia), 'NCIT:C3224' (Melanoma), 'NCIT:C4872' (Lung Adenocarcinoma), 'NCIT:C3167' (Neuroblastoma). |
| `limit` | integer | no | Maximum number of biosamples to return (default: 10). Default `10` |
| `skip` | integer | no | Number of results to skip for pagination (default: 0). Default `0` |

??? note "Upstream description"

    Search cancer tumor biosamples in the Progenetix database by NCIt disease ontology code. Progenetix contains genome-wide copy number variation (CNV) profiles from over 100,000 tumor samples. Returns biosample IDs, histological diagnosis, pathological stage, and external references. Uses the GA4GH Beacon v2 protocol. Example: NCIT:C3058 returns 28 glioblastoma biosamples with histological classification and study references.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Progenetix_search_biosamples',
     'arguments': {'filters': '...'}},
)
```

## `SCP_get_study` {#SCP_get_study}

Retrieves one Broad Single Cell Portal study by accession with its full description and metadata.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** single-cell study, broad institute, scrna-seq, study accession, dataset
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | SCP accession, e.g. 'SCP1'. Case-insensitive. |

??? note "Upstream description"

    Retrieve one Single Cell Portal study by accession, with its full untruncated description. Accessions look like 'SCP1' and are returned by SCP_search_studies and SCP_list_studies. Private or detached studies are not exposed by the public API and will return an error.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SCP_get_study', 'arguments': {'accession': '...'}},
)
```

## `SCP_list_studies` {#SCP_list_studies}

Lists public Single Cell Portal studies sorted by cell count, optionally filtered by minimum cells or keyword.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** single-cell studies, dataset catalog, cell count, scrna-seq, broad institute
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `keyword` | string \| null | no | Optional case-insensitive filter matched against study title and description, e.g. 'brain'. |
| `min_cells` | integer \| null | no | Only return studies with at least this many cells, e.g. 100000. |
| `limit` | integer \| null | no | Maximum studies to return (default 20, max 100). |
| `description_chars` | integer \| null | no | Truncate each description to this many characters (default 300). Set 0 to omit descriptions. |

??? note "Upstream description"

    List the public Single Cell Portal study catalog (1000+ studies), sorted by cell count descending. Optionally filter by a minimum cell count or a keyword matched against title and description. Use this to find the largest available single-cell studies, or to survey what the portal holds. Returns accession, title, cell count, gene count, and study URL.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SCP_list_studies', 'arguments': {}},
)
```

## `SCP_search_studies` {#SCP_search_studies}

Searches Broad Single Cell Portal studies by keyword and returns accession, title, cell and gene counts and URL.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `query`
- **Keywords:** single-cell studies, dataset search, scrna-seq, broad institute, cell count
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | yes | Search terms, e.g. 'lung', 'glioblastoma', 'COVID-19', 'pancreatic islet'. Case-insensitive. |
| `limit` | integer \| null | no | Maximum studies to return (default 20, max 100). |
| `page` | integer \| null | no | 1-based result page. Use with the total_pages value returned in metadata. |
| `description_chars` | integer \| null | no | Truncate each description to this many characters (default 500). Set 0 to omit descriptions. |

??? note "Upstream description"

    Keyword search across the Broad Institute Single Cell Portal, which hosts 1000+ single-cell studies covering 80M+ cells. Searches study titles, descriptions, and metadata. Returns accession, title, cell count, gene count, and study URL. Example: query='glioblastoma'. Complements CxGDisc_search_datasets (CELLxGENE) and SCXA_list_experiments (EBI), which index different study sets. Descriptions are truncated; use SCP_get_study for the full text.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SCP_search_studies', 'arguments': {'query': '...'}},
)
```

## `SCXA_list_experiments` {#SCXA_list_experiments}

Lists curated single-cell RNA-seq experiments in the EBI Single Cell Expression Atlas with species, cell counts, technology and experimental factors.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** scrna-seq datasets, single-cell experiments, 10x, smart-seq2, human cell atlas, expression atlas
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `species` | string | no | Species name to filter by (e.g., 'Homo sapiens', 'Mus musculus', 'Drosophila melanogaster'). Case-insensitive. |
| `keyword` | string | no | Keyword to search in experiment descriptions and factors (e.g., 'pancreas', 'brain', 'cancer', 'immune'). |
| `limit` | integer | no | Maximum number of experiments to return (default 20, max 100). Default `20` |

??? note "Upstream description"

    List single-cell RNA-seq experiments from EBI Single Cell Expression Atlas. Returns 380+ curated scRNA-seq datasets with species, cell counts, technology (10x, Smart-seq2, Drop-seq), and experimental factors (cell type, tissue, disease). Filter by species or keyword. Covers Human Cell Atlas, GTEx, Tabula Sapiens, and more. Example: filtering by 'Homo sapiens' returns 159 experiments.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SCXA_list_experiments', 'arguments': {}},
)
```

## `SRA_get_experiment` {#SRA_get_experiment}

Gets metadata for an SRA sequencing experiment by accession, including organism, platform, library details, study and run statistics.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `accession`
- **Keywords:** sequence read archive, srx, library strategy, sequencing platform, runs
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `accession` | string | yes | SRA accession number (e.g., 'SRX26927637', 'SRP000001', 'SRS000001') |

??? note "Upstream description"

    Get detailed metadata for a specific SRA experiment by accession number. Returns experiment title, organism, sequencing platform, library information, study details, and run statistics. Accepts SRA experiment accessions (SRX/ERX/DRX), study accessions (SRP/ERP/DRP), or sample accessions (SRS/ERS/DRS).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SRA_get_experiment', 'arguments': {'accession': '...'}},
)
```

## `SRA_search_experiments` {#SRA_search_experiments}

Searches NCBI SRA for sequencing experiments by keyword, organism, library strategy or platform.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** sequence read archive, raw reads, illumina, library strategy, rna-seq
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `query` | string | no | Search keywords (e.g., 'human gut metagenome', 'BRCA1 RNA-Seq', 'single cell lung cancer') |
| `organism` | string | no | Organism name to filter by (e.g., 'Homo sapiens', 'Mus musculus', 'human gut metagenome') |
| `library_strategy` | string | no | Sequencing strategy filter (e.g., 'RNA-Seq', 'WGS', 'WXS', 'ChIP-Seq', 'ATAC-seq', 'Bisulfite-Seq', 'AMPLICON') |
| `platform` | string | no | Sequencing platform filter (e.g., 'ILLUMINA', 'PACBIO_SMRT', 'OXFORD_NANOPORE', 'ION_TORRENT') |
| `limit` | integer | no | Maximum results to return (1-50, default 10) Default `10` |

??? note "Upstream description"

    Search the NCBI Sequence Read Archive (SRA) for sequencing experiments by keyword, organism, library strategy, or sequencing platform. SRA is the largest repository of raw sequencing data, containing data from Illumina, PacBio, Oxford Nanopore, and other platforms. Library strategies include: RNA-Seq, WGS (whole genome), WXS (whole exome), ChIP-Seq, ATAC-seq, Bisulfite-Seq, AMPLICON, scRNA-Seq, and more.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SRA_search_experiments', 'arguments': {}},
)
```

## `UCSCCellBrowser_get_dataset` {#UCSCCellBrowser_get_dataset}

Gets details of a UCSC Cell Browser single-cell dataset including abstract, organisms, body parts, diseases and child datasets.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `name`
- **Keywords:** single-cell dataset, cell browser, scrna-seq, dataset metadata
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `name` | string | yes | Dataset identifier, e.g. 'organoid-22q11'. Returned by UCSCCellBrowser_search_datasets. |

??? note "Upstream description"

    Get detail for one UCSC Cell Browser dataset by name, including its abstract, organisms, body parts, diseases, parent collections, and any child datasets. Dataset names look like 'organoid-22q11' and come from UCSCCellBrowser_search_datasets. For collection entries, the child_datasets array lists the individual datasets inside.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UCSCCellBrowser_get_dataset',
     'arguments': {'name': '...'}},
)
```

## `UCSCCellBrowser_list_facets` {#UCSCCellBrowser_list_facets}

Lists distinct organism, body part or disease facet values with dataset counts in the UCSC Cell Browser catalog.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** `facet`
- **Keywords:** facets, organisms, body parts, single-cell catalog
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `facet` | string | yes | Which facet to enumerate. One of `organisms`, `body_parts`, `diseases` |

??? note "Upstream description"

    List the distinct values available for a UCSC Cell Browser filter facet, with the number of datasets carrying each value. Valid facets are 'organisms', 'body_parts', and 'diseases'. Call this before UCSCCellBrowser_search_datasets to discover exact filter values, since the catalog uses controlled vocabulary like 'Human (H. sapiens)'.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UCSCCellBrowser_list_facets',
     'arguments': {'facet': 'organisms'}},
)
```

## `UCSCCellBrowser_search_datasets` {#UCSCCellBrowser_search_datasets}

Searches the UCSC Cell Browser catalog of single-cell datasets by organism, body part, disease or keyword.

- **Category:** `sequencing_data` (Sequencing datasets, runs, experiments and read files)
- **Required arguments:** none
- **Keywords:** single-cell datasets, scrna-seq, cell browser, dataset search
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `organism` | string \| null | no | Substring match on organism, e.g. 'Human', 'Mouse', 'Zebrafish'. |
| `body_part` | string \| null | no | Substring match on body part, e.g. 'brain', 'lung', 'retina', 'blood'. |
| `disease` | string \| null | no | Substring match on disease, e.g. 'cancer', 'COVID-19', 'Healthy'. |
| `keyword` | string \| null | no | Substring match on dataset name or display label, e.g. 'organoid'. |
| `limit` | integer \| null | no | Maximum datasets to return (default 25, max 100). |

??? note "Upstream description"

    Search the UCSC Cell Browser catalog of 300+ curated single-cell datasets, filtering by organism, body part, disease, or keyword. Each dataset is consistently annotated, so filters are reliable across studies. Example: organism='Human', body_part='brain'. Returns dataset name, label, organisms, body parts, diseases, and browser URL. Use UCSCCellBrowser_list_facets first to see valid filter values.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'UCSCCellBrowser_search_datasets', 'arguments': {}},
)
```
