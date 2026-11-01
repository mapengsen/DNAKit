# 10x Genomics Cloud MCP (`10x-genomics`)

DNAKit's MCP tool index keeps **3** DNA/gene tools of this server (27 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [10XGenomics/txg-mcp @ `2013be5`](https://github.com/10XGenomics/txg-mcp/tree/2013be5d9c488f7914462cae0cf873939f6dd6b9) |
| License | MIT (MCP connector); separate 10x Cloud and TXG CLI terms apply |
| Command | `node server/index.js` |
| Transports | stdio |
| Installation | Separately prepare and unpack the upstream 10x Genomics MCP bundle, including its server/, node_modules/ and platform bin/ files. Set workdir to that unpacked bundle, executable_path to installed Node.js &gt;=18, and the named TXG_CLI_ACCESS_TOKEN environment variable. DNAKit does not download a bundle or TXG binary. |
| Required environment | `TXG_CLI_ACCESS_TOKEN` (secret) |
| Effects | Selected tools can upload sequencing files, read or modify cloud projects and references, and start paid analyses. Upstream only sets accept-payment=true when that tool argument is explicitly supplied. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`create_cellranger_aggr_analysis`](#create_cellranger_aggr_analysis) | `expression` | Creates a 10x Cloud Cell Ranger aggr analysis that aggregates multiple count/multi runs from a config CSV, with optional depth normalization. |
| [`create_cellranger_count_analysis`](#create_cellranger_count_analysis) | `expression` | Creates a 10x Cloud Cell Ranger count analysis that aligns single-sample FASTQs to a reference transcriptome and counts gene expression per cell. |
| [`create_cellranger_multi_analysis`](#create_cellranger_multi_analysis) | `expression` | Creates a 10x Cloud Cell Ranger multi analysis for Flex, multiplexed, feature-barcoding or multimodal experiments from a config CSV. |

## `create_cellranger_aggr_analysis` {#create_cellranger_aggr_analysis}

Creates a 10x Cloud Cell Ranger aggr analysis that aggregates multiple count/multi runs from a config CSV, with optional depth normalization.

- **Category:** `expression` (Gene-level expression and transcript abundance)
- **Required arguments:** `analysis_name`, `csv_path`
- **Keywords:** cell ranger aggr, single-cell rna-seq, gene-barcode matrix, aggregation, 10x cloud

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `analysis_name` | string | yes | Name of the aggregation analysis |
| `csv_path` | string | yes | Path to aggr config CSV file. CRITICAL: if you need to generate the file, you must use 'get_aggr_csv_config_spec' tool to get the correct format specification and examples. |
| `project_id` | string \| null | no | ID of existing project to create the analysis in Default `None` |
| `project_name` | string \| null | no | Name of a new project to create for this analysis (alternative to project_id) Default `None` |
| `normalize` | string \| null | no | Normalization method: 'mapped' (default) or 'none' Default `mapped` |
| `description` | string \| null | no | Description for the aggregation analysis Default `None` |
| `accept_payment` | string \| null | no | CRITICAL: Always omit this parameter on first attempt as most analyses do not require payment. ONLY set to "true" when: (1) analysis creation failed with payment error, (2) you showed the user the cost, (3) user explicitly confirmed. Never preemptively ask about payment. Default `None` |

??? note "Upstream description"

    Creates a new Cell Ranger 'aggr' analysis to aggregate multiple runs. CRITICAL: Always explicitly confirm parameters before creating analysis (can be expensive, long-running compute). Ask about experimental setup, explain parameter choices, show final parameters in table, and ASK 'Please review and confirm these parameters before I create the analysis.'

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    '10x-genomics',
    'create_cellranger_aggr_analysis',
    {'analysis_name': '...', 'csv_path': '...'},
)
```

## `create_cellranger_count_analysis` {#create_cellranger_count_analysis}

Creates a 10x Cloud Cell Ranger count analysis that aligns single-sample FASTQs to a reference transcriptome and counts gene expression per cell.

- **Category:** `expression` (Gene-level expression and transcript abundance)
- **Required arguments:** `analysis_name`, `transcriptome`, `fastqs`
- **Keywords:** cell ranger count, single-cell rna-seq, fastq, transcriptome alignment, gene-barcode matrix, 10x cloud

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `analysis_name` | string | yes | Name of the analysis to create |
| `transcriptome` | string | yes | Reference transcriptome to use (e.g., 'refdata-cellranger-GRCh38-2024-A', or custom reference ID). Use the 'list_custom_references' or 'list_prebuilt_references' tools to find available references. |
| `fastqs` | array of string | yes | List of FASTQ file paths or FASTQ URIs to analyze. Can be local paths or previously uploaded FASTQs in the format txg://fastqs/&lt;filename&gt; or txg://fastqs/&lt;fastq_uuid&gt; |
| `project_id` | string \| null | no | ID of existing project to create the analysis in Default `None` |
| `project_name` | string \| null | no | Name of a new project to create for this analysis (alternative to project_id) Default `None` |
| `expect_cells` | integer \| null | no | Expected number of recovered cells, used as input to cell calling algorithm. Mutually exclusive with force_cells. Default `None` |
| `force_cells` | integer \| null | no | Force pipeline to use this number of cells, bypassing cell calling algorithm. Minimum: 10. Mutually exclusive with expect_cells. Default `None` |
| `chemistry` | string \| null | no | Assay chemistry version (default: 'auto' for automatic detection) Default `auto` |
| `cell_annotation_model` | string \| null | no | Cell annotation model to use. If not provided, does not run cell annotation. Use 'auto' for default model for the species. or specify a particular model ID from 'list_annotation_models' tool. Default `None` |
| `include_introns` | boolean | no | Whether to include intronic reads in count (default: true) Default `True` |
| `create_bam` | boolean | no | Enable or disable BAM file generation (default: true). Setting to false reduces computation time and output size. Default `True` |
| `nosecondary` | boolean | no | Disable secondary analysis (e.g., clustering) (default: false) Default `False` |
| `r1_length` | integer \| null | no | Hard trim the input Read 1 to this length before analysis Default `None` |
| `r2_length` | integer \| null | no | Hard trim the input Read 2 to this length before analysis Default `None` |
| `accept_payment` | string \| null | no | CRITICAL: Always omit this parameter on first attempt as most analyses do not require payment. ONLY set to "true" when: (1) analysis creation failed with payment error, (2) you showed the user the cost, (3) user explicitly confirmed. Never preemptively ask about payment. Default `None` |

??? note "Upstream description"

    Creates a new Cell Ranger 'count' analysis. ONLY for simple single-sample gene expression without additional modalities. NOT for Flex assays. Pipeline selection: Use 'multi' for Flex assays, multimodal experiments, feature barcoding, or multiplexed samples. Use 'count' for simple single-sample gene expression only. If unsure, ask user about assay type. CRITICAL: Always explicitly confirm parameters before creating analysis (can be expensive, long-running compute). Ask about experimental setup, explain parameter choices, show final parameters in table, and ASK 'Please review and confirm these parameters before I create the analysis.'

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    '10x-genomics',
    'create_cellranger_count_analysis',
    {'analysis_name': '...', 'transcriptome': '...', 'fastqs': []},
)
```

## `create_cellranger_multi_analysis` {#create_cellranger_multi_analysis}

Creates a 10x Cloud Cell Ranger multi analysis for Flex, multiplexed, feature-barcoding or multimodal experiments from a config CSV.

- **Category:** `expression` (Gene-level expression and transcript abundance)
- **Required arguments:** `analysis_name`, `csv_path`
- **Keywords:** cell ranger multi, flex, feature barcoding, multiplexing, single-cell rna-seq, 10x cloud

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `analysis_name` | string | yes | Name of the analysis to create |
| `csv_path` | string | yes | Path to the multi config CSV file. CRITICAL: if you need to generate the file, you must use 'get_multi_csv_config_spec' tool to get the correct format specification and examples. |
| `project_id` | string \| null | no | ID of existing project to create the analysis in. Either project_id or project_name must be specified Default `None` |
| `project_name` | string \| null | no | Name of a new project to create for this analysis (alternative to project_id) Default `None` |
| `description` | string \| null | no | Analysis description Default `None` |
| `accept_payment` | string \| null | no | CRITICAL: Always omit this parameter on first attempt as most analyses do not require payment. ONLY set to "true" when: (1) analysis creation failed with payment error, (2) you showed the user the cost, (3) user explicitly confirmed. Never preemptively ask about payment. Default `None` |

??? note "Upstream description"

    Creates a new Cell Ranger 'multi' analysis. REQUIRED for Flex assays, multimodal experiments, feature barcoding, multiplexed samples. IMPORTANT: if csv file is not already provided, use 'get_multi_csv_config_spec' tool FIRST to get the correct CSV format for 10x Cloud. Pipeline selection: Use 'multi' for Flex assays, multimodal experiments, feature barcoding, or multiplexed samples. Use 'count' for simple single-sample gene expression only. If unsure, ask user about assay type. CRITICAL: Always explicitly confirm parameters before creating analysis (can be expensive, long-running compute). Ask about experimental setup, explain parameter choices, show final parameters in table, and ASK 'Please review and confirm these parameters before I create the analysis.'

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    '10x-genomics',
    'create_cellranger_multi_analysis',
    {'analysis_name': '...', 'csv_path': '...'},
)
```
