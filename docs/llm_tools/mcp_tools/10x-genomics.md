# 10x Genomics Cloud MCP (`10x-genomics`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **3** 个 DNA/基因相关工具（服务器共列出 27 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [10XGenomics/txg-mcp @ `2013be5`](https://github.com/10XGenomics/txg-mcp/tree/2013be5d9c488f7914462cae0cf873939f6dd6b9) |
| 许可证 | MIT (MCP connector); separate 10x Cloud and TXG CLI terms apply |
| 启动命令 | `node server/index.js` |
| 传输方式 | stdio |
| 安装 | Separately prepare and unpack the upstream 10x Genomics MCP bundle, including its server/, node_modules/ and platform bin/ files. Set workdir to that unpacked bundle, executable_path to installed Node.js &gt;=18, and the named TXG_CLI_ACCESS_TOKEN environment variable. DNAKit does not download a bundle or TXG binary. |
| 必需环境变量 | `TXG_CLI_ACCESS_TOKEN`（保密） |
| 副作用 | Selected tools can upload sequencing files, read or modify cloud projects and references, and start paid analyses. Upstream only sets accept-payment=true when that tool argument is explicitly supplied. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`create_cellranger_aggr_analysis`](#create_cellranger_aggr_analysis) | `expression` | Creates a 10x Cloud Cell Ranger aggr analysis that aggregates multiple count/multi runs from a config CSV, with optional depth normalization. |
| [`create_cellranger_count_analysis`](#create_cellranger_count_analysis) | `expression` | Creates a 10x Cloud Cell Ranger count analysis that aligns single-sample FASTQs to a reference transcriptome and counts gene expression per cell. |
| [`create_cellranger_multi_analysis`](#create_cellranger_multi_analysis) | `expression` | Creates a 10x Cloud Cell Ranger multi analysis for Flex, multiplexed, feature-barcoding or multimodal experiments from a config CSV. |

## `create_cellranger_aggr_analysis` {#create_cellranger_aggr_analysis}

Creates a 10x Cloud Cell Ranger aggr analysis that aggregates multiple count/multi runs from a config CSV, with optional depth normalization.

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`analysis_name`, `csv_path`
- **关键词：**cell ranger aggr, single-cell rna-seq, gene-barcode matrix, aggregation, 10x cloud

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `analysis_name` | string | 是 | Name of the aggregation analysis |
| `csv_path` | string | 是 | Path to aggr config CSV file. CRITICAL: if you need to generate the file, you must use 'get_aggr_csv_config_spec' tool to get the correct format specification and examples. |
| `project_id` | string \| null | 否 | ID of existing project to create the analysis in 默认 `None` |
| `project_name` | string \| null | 否 | Name of a new project to create for this analysis (alternative to project_id) 默认 `None` |
| `normalize` | string \| null | 否 | Normalization method: 'mapped' (default) or 'none' 默认 `mapped` |
| `description` | string \| null | 否 | Description for the aggregation analysis 默认 `None` |
| `accept_payment` | string \| null | 否 | CRITICAL: Always omit this parameter on first attempt as most analyses do not require payment. ONLY set to "true" when: (1) analysis creation failed with payment error, (2) you showed the user the cost, (3) user explicitly confirmed. Never preemptively ask about payment. 默认 `None` |

??? note "上游说明（英文原文）"

    Creates a new Cell Ranger 'aggr' analysis to aggregate multiple runs. CRITICAL: Always explicitly confirm parameters before creating analysis (can be expensive, long-running compute). Ask about experimental setup, explain parameter choices, show final parameters in table, and ASK 'Please review and confirm these parameters before I create the analysis.'

安装服务器后通过 DNAKit 调用：

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

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`analysis_name`, `transcriptome`, `fastqs`
- **关键词：**cell ranger count, single-cell rna-seq, fastq, transcriptome alignment, gene-barcode matrix, 10x cloud

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `analysis_name` | string | 是 | Name of the analysis to create |
| `transcriptome` | string | 是 | Reference transcriptome to use (e.g., 'refdata-cellranger-GRCh38-2024-A', or custom reference ID). Use the 'list_custom_references' or 'list_prebuilt_references' tools to find available references. |
| `fastqs` | array of string | 是 | List of FASTQ file paths or FASTQ URIs to analyze. Can be local paths or previously uploaded FASTQs in the format txg://fastqs/&lt;filename&gt; or txg://fastqs/&lt;fastq_uuid&gt; |
| `project_id` | string \| null | 否 | ID of existing project to create the analysis in 默认 `None` |
| `project_name` | string \| null | 否 | Name of a new project to create for this analysis (alternative to project_id) 默认 `None` |
| `expect_cells` | integer \| null | 否 | Expected number of recovered cells, used as input to cell calling algorithm. Mutually exclusive with force_cells. 默认 `None` |
| `force_cells` | integer \| null | 否 | Force pipeline to use this number of cells, bypassing cell calling algorithm. Minimum: 10. Mutually exclusive with expect_cells. 默认 `None` |
| `chemistry` | string \| null | 否 | Assay chemistry version (default: 'auto' for automatic detection) 默认 `auto` |
| `cell_annotation_model` | string \| null | 否 | Cell annotation model to use. If not provided, does not run cell annotation. Use 'auto' for default model for the species. or specify a particular model ID from 'list_annotation_models' tool. 默认 `None` |
| `include_introns` | boolean | 否 | Whether to include intronic reads in count (default: true) 默认 `True` |
| `create_bam` | boolean | 否 | Enable or disable BAM file generation (default: true). Setting to false reduces computation time and output size. 默认 `True` |
| `nosecondary` | boolean | 否 | Disable secondary analysis (e.g., clustering) (default: false) 默认 `False` |
| `r1_length` | integer \| null | 否 | Hard trim the input Read 1 to this length before analysis 默认 `None` |
| `r2_length` | integer \| null | 否 | Hard trim the input Read 2 to this length before analysis 默认 `None` |
| `accept_payment` | string \| null | 否 | CRITICAL: Always omit this parameter on first attempt as most analyses do not require payment. ONLY set to "true" when: (1) analysis creation failed with payment error, (2) you showed the user the cost, (3) user explicitly confirmed. Never preemptively ask about payment. 默认 `None` |

??? note "上游说明（英文原文）"

    Creates a new Cell Ranger 'count' analysis. ONLY for simple single-sample gene expression without additional modalities. NOT for Flex assays. Pipeline selection: Use 'multi' for Flex assays, multimodal experiments, feature barcoding, or multiplexed samples. Use 'count' for simple single-sample gene expression only. If unsure, ask user about assay type. CRITICAL: Always explicitly confirm parameters before creating analysis (can be expensive, long-running compute). Ask about experimental setup, explain parameter choices, show final parameters in table, and ASK 'Please review and confirm these parameters before I create the analysis.'

安装服务器后通过 DNAKit 调用：

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

- **类别：**`expression`（Gene-level expression and transcript abundance）
- **必填参数：**`analysis_name`, `csv_path`
- **关键词：**cell ranger multi, flex, feature barcoding, multiplexing, single-cell rna-seq, 10x cloud

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `analysis_name` | string | 是 | Name of the analysis to create |
| `csv_path` | string | 是 | Path to the multi config CSV file. CRITICAL: if you need to generate the file, you must use 'get_multi_csv_config_spec' tool to get the correct format specification and examples. |
| `project_id` | string \| null | 否 | ID of existing project to create the analysis in. Either project_id or project_name must be specified 默认 `None` |
| `project_name` | string \| null | 否 | Name of a new project to create for this analysis (alternative to project_id) 默认 `None` |
| `description` | string \| null | 否 | Analysis description 默认 `None` |
| `accept_payment` | string \| null | 否 | CRITICAL: Always omit this parameter on first attempt as most analyses do not require payment. ONLY set to "true" when: (1) analysis creation failed with payment error, (2) you showed the user the cost, (3) user explicitly confirmed. Never preemptively ask about payment. 默认 `None` |

??? note "上游说明（英文原文）"

    Creates a new Cell Ranger 'multi' analysis. REQUIRED for Flex assays, multimodal experiments, feature barcoding, multiplexed samples. IMPORTANT: if csv file is not already provided, use 'get_multi_csv_config_spec' tool FIRST to get the correct CSV format for 10x Cloud. Pipeline selection: Use 'multi' for Flex assays, multimodal experiments, feature barcoding, or multiplexed samples. Use 'count' for simple single-sample gene expression only. If unsure, ask user about assay type. CRITICAL: Always explicitly confirm parameters before creating analysis (can be expensive, long-running compute). Ask about experimental setup, explain parameter choices, show final parameters in table, and ASK 'Please review and confirm these parameters before I create the analysis.'

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    '10x-genomics',
    'create_cellranger_multi_analysis',
    {'analysis_name': '...', 'csv_path': '...'},
)
```
