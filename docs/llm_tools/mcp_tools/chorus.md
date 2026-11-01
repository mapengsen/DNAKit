# Chorus (`chorus`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **17** 个 DNA/基因相关工具（服务器共列出 24 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [pinellolab/chorus @ `cc8aabe`](https://github.com/pinellolab/chorus/tree/cc8aabe9f5e4e562f0fc8b8bfa838fb100cba52b) |
| 许可证 | MIT |
| 启动命令 | `chorus-mcp` |
| 传输方式 | stdio |
| 安装 | Install the pinellolab/chorus source project and provision its oracle environments, model weights and references separately. The PyPI package named chorus is unrelated. Set executable_path to the installed chorus-mcp. |
| 副作用 | Predictions can consume GPU/CPU, write reports, and download missing model assets inside the upstream server. Pre-provision assets when offline. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`analyze_region_swap`](#analyze_region_swap) | `prediction_model` | Replaces a genomic region with a custom DNA sequence and scores wild-type vs replacement predictions across chromatin, TF, histone, CAGE and RNA layers. |
| [`analyze_variant_multilayer`](#analyze_variant_multilayer) | `prediction_model` | Analyzes a variant's predicted regulatory impact across molecular layers with modality-specific scoring and nearby-gene expression effects. |
| [`discover_variant`](#discover_variant) | `prediction_model` | Predicts a variant's effect across all oracle tracks, ranks tracks and cell types by effect, and returns a report. |
| [`discover_variant_cell_types`](#discover_variant_cell_types) | `prediction_model` | Screens DNASE/ATAC tracks to rank cell types affected by a variant, then runs multi-layer analysis for the top cell types. |
| [`fine_map_causal_variant`](#fine_map_causal_variant) | `prediction_model` | Scores a sentinel GWAS variant and its LD proxies across regulatory layers and ranks them by a composite causal score. |
| [`get_gene_tss`](#get_gene_tss) | `genome_annotation` | Gets transcription start site (TSS) positions for a gene symbol from the reference annotation. |
| [`get_genes_in_region`](#get_genes_in_region) | `genome_annotation` | Gets gene annotations overlapping a chromosome region. |
| [`predict`](#predict) | `prediction_model` | Makes a wild-type prediction of selected assay tracks for a genomic region, returning summary stats and bedgraph files. |
| [`predict_region_insertion`](#predict_region_insertion) | `prediction_model` | Inserts a DNA sequence at a genomic position and predicts assay track activity with the loaded oracle. |
| [`predict_region_replacement`](#predict_region_replacement) | `prediction_model` | Replaces a genomic region with a custom DNA sequence and predicts assay track activity. |
| [`predict_variant_effect`](#predict_variant_effect) | `prediction_model` | Predicts a variant's effect by comparing reference and alternate allele predictions across selected or all tracks. |
| [`predict_variant_effect_on_gene`](#predict_variant_effect_on_gene) | `prediction_model` | Predicts a variant's effect on a nearby gene's expression using CAGE TSS and RNA exon quantification. |
| [`score_ism`](#score_ism) | `prediction_model` | Performs in-silico saturation mutagenesis in a window around a position, returning per-base importance for motif logos. |
| [`score_prediction_region`](#score_prediction_region) | `prediction_model` | Predicts a region and summarises signal within a specified sub-region using mean, max, sum or median. |
| [`score_variant_batch`](#score_variant_batch) | `prediction_model` | Scores multiple variants through multi-layer analysis and returns a table ranked by effect magnitude. |
| [`score_variant_effect_at_region`](#score_variant_effect_at_region) | `prediction_model` | Predicts a variant effect and scores ref/alt signal in a chosen sub-region or a window around the variant. |
| [`simulate_integration`](#simulate_integration) | `prediction_model` | Simulates inserting a construct sequence at a genomic position and scores disruption across regulatory layers. |

## `analyze_region_swap` {#analyze_region_swap}

Replaces a genomic region with a custom DNA sequence and scores wild-type vs replacement predictions across chromatin, TF, histone, CAGE and RNA layers.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `region`, `replacement_sequence`, `assay_ids`
- **关键词：**promoter swap, enhancer replacement, sequence-to-function, alphagenome, enformer

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `region` | string | 是 | Region to replace as "chr1:1000000-1001000". |
| `replacement_sequence` | string | 是 | DNA sequence to insert in place of the region. |
| `assay_ids` | array of string | 是 | List of assay identifiers for multi-layer scoring. |
| `gene_name` | string \| null | 否 | Optional gene for expression scoring. 默认 `None` |
| `description` | string \| null | 否 | Optional short label for the swap (e.g. "SV40 promoter"). 默认 `None` |
| `user_prompt` | string \| null | 否 | Original user prompt, rendered at the top of the report. 默认 `None` |

??? note "上游说明（英文原文）"

    Replace a genomic region with a custom sequence and score effects across all layers.

    Compares wild-type vs replacement predictions using the same multi-layer
    scoring as variant analysis (chromatin, TF binding, histone, CAGE, RNA).

    Use cases: promoter swaps, enhancer replacements, regulatory element engineering.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'analyze_region_swap',
    {'oracle_name': '...',
     'region': '...',
     'replacement_sequence': '...',
     'assay_ids': []},
)
```

## `analyze_variant_multilayer` {#analyze_variant_multilayer}

Analyzes a variant's predicted regulatory impact across molecular layers with modality-specific scoring and nearby-gene expression effects.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `position`, `ref_allele`, `alt_alleles`, `assay_ids`
- **关键词：**variant effect prediction, chromatin, tf binding, cage, log2 fold change

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `position` | string | 是 | Variant position as ``"chr1:1050000"`` OR an rsID like ``"rs12740374"`` (requires LDlink token; resolves to coords via the LDproxy API). |
| `ref_allele` | string | 是 | Reference allele. |
| `alt_alleles` | array of string | 是 | Alternate alleles. |
| `assay_ids` | array of string | 是 | List of assay identifiers covering different layers (e.g. DNASE, CAGE, ChIP tracks for multi-layer coverage). Pass an empty list or None to score all tracks on oracles that support it (AlphaGenome, Enformer, Borzoi). |
| `gene_name` | string \| null | 否 | Gene symbol for RNA expression scoring (e.g. "SORT1"). If omitted, the nearest gene is auto-detected. 默认 `None` |
| `region` | string \| null | 否 | Input region as "chr1:1000000-1393216". If omitted, auto-centered. 默认 `None` |
| `igv_raw` | boolean | 否 | When True, the IGV browser in the HTML report shows raw signal with autoscale instead of the layer-aware rescaled view. Table scores are unaffected. 默认 `False` |
| `ldlink_token` | string \| null | 否 | LDlink API token (only used when ``position`` is an rsID). Register free at https://ldlink.nih.gov/?tab=apiaccess or set ``LDLINK_TOKEN``. 默认 `None` |
| `genome_build` | string | 否 | Reference build for the rsID lookup (``"grch38"`` / ``"hg38"`` default; or ``"grch37"`` / ``"hg19"``). Only used when ``position`` is an rsID. 默认 `grch38` |
| `user_prompt` | string \| null | 否 | The user's original natural-language question. Claude should forward this verbatim whenever calling from an MCP conversation — it is rendered at the top of the report for traceability. 默认 `None` |

??? note "上游说明（英文原文）"

    Analyze a variant's regulatory impact across all molecular layers.

    Scores each track using modality-specific strategies:
    - Chromatin (DNASE/ATAC): log2 fold-change of sum in 501bp window
    - TF binding (ChIP-TF): log2 fold-change of sum in 501bp window
    - Histone marks (ChIP-Histone): log2 fold-change of sum in 2001bp window
    - TSS activity (CAGE): log2 fold-change of sum in 501bp window
    - Gene expression (RNA): log fold-change of mean over gene exons
    - Promoter activity (MPRA): simple difference

    For non-coding variants, nearby genes are auto-detected within the
    prediction window so that RNA expression effects can be scored even
    without an explicit gene_name.

    Returns a structured report with scores organized by regulatory layer,
    plus a markdown summary for interpretation. Every report carries the
    original user prompt at the top so it stays interpretable months later.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'analyze_variant_multilayer',
    {'oracle_name': '...',
     'position': '...',
     'ref_allele': '...',
     'alt_alleles': [],
     'assay_ids': []},
)
```

## `discover_variant` {#discover_variant}

Predicts a variant's effect across all oracle tracks, ranks tracks and cell types by effect, and returns a report.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `position`, `ref_allele`, `alt_alleles`
- **关键词：**variant interpretation, track ranking, alphagenome, borzoi, cell types

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `position` | string | 是 | Variant position as ``"chr1:109274968"`` OR an rsID like ``"rs12740374"`` (requires LDlink token; resolves to coords via the LDproxy API). |
| `ref_allele` | string | 是 | Reference allele (e.g. "G"). |
| `alt_alleles` | array of string | 是 | Alternate alleles (e.g. \["T"\]). |
| `gene_name` | string \| null | 否 | Optional gene for expression analysis. 默认 `None` |
| `top_n` | integer | 否 | Number of top tracks per regulatory layer to show. 默认 `3` |
| `igv_raw` | boolean | 否 | 默认 `False` |
| `user_prompt` | string \| null | 否 | Original user prompt, forwarded into the report header. 默认 `None` |
| `ranking_metric` | string | 否 | How to rank tracks and cell types. Default ``"alt_x_abs_effect"`` (``alt_value × \|log2FC\|``) avoids inflating closed-baseline tracks where a small absolute change yields a large fold-change. Use ``"abs_effect"`` for the historical raw \|log2FC\| ranking. 默认 `alt_x_abs_effect` |
| `min_ref_value` | number | 否 | Threshold used by ``"abs_effect_min_ref"``. 默认 `0.0` |
| `ldlink_token` | string \| null | 否 | LDlink API token (only used when ``position`` is an rsID). Register free at https://ldlink.nih.gov/?tab=apiaccess. 默认 `None` |
| `genome_build` | string | 否 | Reference build for the rsID lookup (``"grch38"`` / ``"hg38"`` default; or ``"grch37"`` / ``"hg19"``). 默认 `grch38` |

??? note "上游说明（英文原文）"

    Discover which cell types and regulatory layers are most affected by a variant.

    Predicts variant effect across ALL available tracks (thousands for
    Enformer/Borzoi/AlphaGenome, or iterates all models for ChromBPNet/LegNet),
    ranks by ``ranking_metric``, and returns the top hits with a full
    report.

    This is the primary tool for variant interpretation — it tells you WHERE
    the variant has impact without requiring you to pre-select tracks.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'discover_variant',
    {'oracle_name': '...',
     'position': '...',
     'ref_allele': '...',
     'alt_alleles': []},
)
```

## `discover_variant_cell_types` {#discover_variant_cell_types}

Screens DNASE/ATAC tracks to rank cell types affected by a variant, then runs multi-layer analysis for the top cell types.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `position`, `ref_allele`, `alt_alleles`
- **关键词：**cell-type discovery, dnase, atac, variant effect, alphagenome

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name (ideally AlphaGenome for broadest cell-type coverage). |
| `position` | string | 是 | Variant position as ``"chr1:1050000"`` OR an rsID like ``"rs12740374"`` (requires LDlink token; resolves to coords via the LDproxy API). |
| `ref_allele` | string | 是 | Reference allele. |
| `alt_alleles` | array of string | 是 | Alternate alleles. |
| `gene_name` | string \| null | 否 | Optional gene to focus expression analysis on. 默认 `None` |
| `top_n` | integer | 否 | Number of top cell types to analyze in detail (default 5). 默认 `5` |
| `min_effect` | number | 否 | Minimum \|log2FC\| in DNASE/ATAC to consider a cell type hit (default 0.15). 默认 `0.15` |
| `user_prompt` | string \| null | 否 | Original user prompt, forwarded into each sub-report. 默认 `None` |
| `ranking_metric` | string | 否 | How to rank cell types. Default ``"alt_x_abs_effect"`` (recommended) ranks by ``alt_signal × \|log2FC\|``, which surfaces the cell type where the post-variant element is most active. Use ``"abs_effect"`` to reproduce the historical \|log2FC\|-only ranking (which has a bias toward closed-baseline cell types — see :func:`chorus.analysis.discovery.discover_cell_types`). Use ``"abs_effect_min_ref"`` with… 默认 `alt_x_abs_effect` |
| `min_ref_value` | number | 否 | Threshold for ``"abs_effect_min_ref"`` (ignored otherwise). Default 0 = no filter. 默认 `0.0` |
| `ldlink_token` | string \| null | 否 | LDlink API token (only used when ``position`` is an rsID). 默认 `None` |
| `genome_build` | string | 否 | Reference build for the rsID lookup (``"grch38"`` / ``"hg38"`` default; or ``"grch37"`` / ``"hg19"``). 默认 `grch38` |

??? note "上游说明（英文原文）"

    Discovery mode: find which cell types are most affected by a variant.

    Use this when you don't know which cell type is relevant — let the model
    tell you where the variant matters most.

    **Two-stage analysis:**
    1. Screens all DNASE/ATAC tracks (~472 cell types on AlphaGenome, ~638
    on Enformer) to rank cell types by ``ranking_metric``.
    2. For each top cell type, runs full multi-layer analysis (chromatin,
    TF, histone, CAGE, RNA) limited to that cell type's tracks.

    **Runtime expectations** (AlphaGenome, single A100):
    - Stage 1 screen: ~30–60 s
    - Stage 2 per-cell-type analysis: ~30 s × ``top_n``
    - Typical end-to-end with default ``top_n=5``: 3–4 minutes

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'discover_variant_cell_types',
    {'oracle_name': '...',
     'position': '...',
     'ref_allele': '...',
     'alt_alleles': []},
)
```

## `fine_map_causal_variant` {#fine_map_causal_variant}

Scores a sentinel GWAS variant and its LD proxies across regulatory layers and ranks them by a composite causal score.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `lead_variant`
- **关键词：**fine-mapping, gwas, ld proxies, ldlink, causal variant

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `lead_variant` | string | 是 | Sentinel variant as "rs12740374" or "chr1:109274968 G&gt;T". |
| `ld_variants` | array \| null | 否 | Optional list of LD variant dicts (schema above). If omitted, auto-fetched from LDlink. 默认 `None` |
| `assay_ids` | array \| null | 否 | Track identifiers. Pass an empty list / None to score all tracks (recommended for AlphaGenome). 默认 `None` |
| `gene_name` | string \| null | 否 | Target gene for expression scoring. 默认 `None` |
| `population` | string | 否 | 1000 Genomes population for LD lookup (default CEU). 默认 `CEU` |
| `r2_threshold` | number | 否 | Minimum r² for LD variants (default 0.8). 默认 `0.8` |
| `ldlink_token` | string \| null | 否 | LDlink API token. Register free at https://ldlink.nih.gov/?tab=apiaccess 默认 `None` |
| `ldlink_timeout` | number | 否 | HTTP timeout in seconds for the LDlink request (default 30). Increase for slow networks or large LD blocks. 默认 `30.0` |
| `genome_build` | string | 否 | Reference build for the LDlink LD lookup. Accepts ``"grch38"`` / ``"hg38"`` (default) or ``"grch37"`` / ``"hg19"``. 默认 `grch38` |
| `snvs_only` | boolean | 否 | When True, restrict scoring to single-nucleotide substitutions (drops insertions, deletions, MNVs, and complex multi-base proxies). Default False — indels and multi-allelic LDlink rows are scored by default as of v0.5.5. Set True to reproduce pre-v0.5.5 behavior. 默认 `False` |
| `user_prompt` | string \| null | 否 | Original user prompt, rendered at the top of the report. 默认 `None` |
| `report_top_n` | integer | 否 | Number of top-ranked variants to render full IGV signal tracks for in the HTML report (default 3). Every variant is scored with the same per-track logic; this only controls how many of the top hits get the (more expensive) interactive browser tracks. Set to 0 to skip the IGV signal tracks entirely. 默认 `3` |

??? note "上游说明（英文原文）"

    Prioritize causal variants from a GWAS locus using multi-layer regulatory evidence.

    Given a sentinel GWAS variant and its LD proxies, scores each variant
    across all regulatory layers and ranks by a **composite causal score**
    combining four components (each in \[0, 1\] after min-max normalization):

    1. ``max_effect``    — largest \|log2FC\| across layers (weight 0.35)
    2. ``n_layers``      — count of layers with effect above threshold (0.25)
    3. ``convergence``   — directional agreement across layers (0.20)
    4. ``ref_activity``  — baseline activity of the variant site (0.20)

    A variant with a *strong effect in many layers, all in the same direction,
    in an already-active region* ends up at the top of the ranking.

    **Two modes:**

    - **Auto-fetch**: provide only ``lead_variant`` + ``ldlink_token`` to
    pull LD proxies from LDlink at the given ``population`` / ``r2_threshold``.
    - **Manual**: provide ``ld_variants`` directly. Each dict needs
    ``chrom``, ``pos``, ``ref``, ``alt``, and optional ``id``, ``r2``.

    **Output (per-variant columns in the ranked table):**

    - ``variant_id`` / rsID (★ marks the sentinel)
    - ``r2`` to the sentinel
    - ``max_effect`` (signed log2FC in the strongest layer)
    - ``n_layers_affected`` — 0 means no layer above threshold
    - ``convergence`` — 1.0 = all effects same sign, 0.0 = split
    - ``composite`` — final ranking score; top row is the most likely
    causal variant

    Use ``result\["rankings"\]\[0\]\["per_layer_scores"\]`` to read off *which*
    layers drove the top candidate's score.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'fine_map_causal_variant',
    {'oracle_name': '...', 'lead_variant': '...'},
)
```

## `get_gene_tss` {#get_gene_tss}

Gets transcription start site (TSS) positions for a gene symbol from the reference annotation.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`gene_name`
- **关键词：**tss, transcription start site, gene annotation, promoter

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `gene_name` | string | 是 | Gene symbol (e.g. "GATA1", "TP53"). |

??? note "上游说明（英文原文）"

    Get transcription start site (TSS) positions for a gene.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'get_gene_tss',
    {'gene_name': '...'},
)
```

## `get_genes_in_region` {#get_genes_in_region}

Gets gene annotations overlapping a chromosome region.

- **类别：**`genome_annotation`（Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation）
- **必填参数：**`chrom`, `start`, `end`
- **关键词：**genomic region, gene annotation, chromosome, interval

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chrom` | string | 是 | Chromosome (e.g. "chr1"). |
| `start` | integer | 是 | Region start position. |
| `end` | integer | 是 | Region end position. |

??? note "上游说明（英文原文）"

    Get gene annotations in a genomic region.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'get_genes_in_region',
    {'chrom': '...', 'start': 0, 'end': 0},
)
```

## `predict` {#predict}

Makes a wild-type prediction of selected assay tracks for a genomic region, returning summary stats and bedgraph files.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `region`, `assay_ids`
- **关键词：**sequence-to-function, bedgraph, dnase, cage, genomic region

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `region` | string | 是 | Genomic region as "chr1:1000000-1393216". |
| `assay_ids` | array of string | 是 | List of assay identifiers (e.g. \["DNASE:K562", "CAGE:K562"\]). |

??? note "上游说明（英文原文）"

    Make a wild-type prediction for a genomic region.

    Returns per-track summary stats (mean/max/min/std). For small predictions
    the full values are included inline; for large ones a downsampled preview
    is returned and full data is saved as bedgraph files.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'predict',
    {'oracle_name': '...', 'region': '...', 'assay_ids': []},
)
```

## `predict_region_insertion` {#predict_region_insertion}

Inserts a DNA sequence at a genomic position and predicts assay track activity with the loaded oracle.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `position`, `sequence`, `assay_ids`
- **关键词：**sequence insertion, in silico, sequence-to-function, regulatory activity

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `position` | string | 是 | Insertion point as "chr1:1050000". |
| `sequence` | string | 是 | DNA sequence to insert. |
| `assay_ids` | array of string | 是 | List of assay identifiers. |

??? note "上游说明（英文原文）"

    Insert a sequence at a genomic position and predict activity.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'predict_region_insertion',
    {'oracle_name': '...',
     'position': '...',
     'sequence': '...',
     'assay_ids': []},
)
```

## `predict_region_replacement` {#predict_region_replacement}

Replaces a genomic region with a custom DNA sequence and predicts assay track activity.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `region`, `replacement_sequence`, `assay_ids`
- **关键词：**region replacement, in silico, sequence-to-function, regulatory activity

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `region` | string | 是 | Genomic region to replace as "chr1:1000000-1001000". |
| `replacement_sequence` | string | 是 | DNA sequence to insert in place of the region. |
| `assay_ids` | array of string | 是 | List of assay identifiers. |

??? note "上游说明（英文原文）"

    Replace a genomic region with a custom sequence and predict activity.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'predict_region_replacement',
    {'oracle_name': '...',
     'region': '...',
     'replacement_sequence': '...',
     'assay_ids': []},
)
```

## `predict_variant_effect` {#predict_variant_effect}

Predicts a variant's effect by comparing reference and alternate allele predictions across selected or all tracks.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `position`, `ref_allele`, `alt_alleles`, `assay_ids`
- **关键词：**variant effect prediction, ref alt, log2 fold change, tracks

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `position` | string | 是 | Variant position as "chr1:1050000". |
| `ref_allele` | string | 是 | Reference allele (e.g. "A"). |
| `alt_alleles` | array of string | 是 | Alternate alleles (e.g. \["G", "T"\]). |
| `assay_ids` | array of string | 是 | List of assay identifiers. |
| `region` | string \| null | 否 | Genomic region as "chr1:1000000-1393216". If omitted, auto-centered on the variant position using the oracle's input window. 默认 `None` |

??? note "上游说明（英文原文）"

    Predict the effect of a genetic variant.

    Compares predictions for reference vs alternate alleles and returns
    per-allele effect sizes with summary statistics.

    One forward pass scores the variant across every track in ``assay_ids``
    (or ALL of the oracle's tracks when ``assay_ids`` is omitted), so request
    the full set in a single call — never loop this tool once per track or per
    cell type. To score many variants, use ``score_variant_batch`` or
    ``fine_map_causal_variant`` (each scores a variant exactly once).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'predict_variant_effect',
    {'oracle_name': '...',
     'position': '...',
     'ref_allele': '...',
     'alt_alleles': [],
     'assay_ids': []},
)
```

## `predict_variant_effect_on_gene` {#predict_variant_effect_on_gene}

Predicts a variant's effect on a nearby gene's expression using CAGE TSS and RNA exon quantification.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `position`, `ref_allele`, `alt_alleles`, `gene_name`, `assay_ids`
- **关键词：**gene expression, cage, rna, variant effect, fold change

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `position` | string | 是 | Variant position as "chr1:1050000". |
| `ref_allele` | string | 是 | Reference allele. |
| `alt_alleles` | array of string | 是 | Alternate alleles. |
| `gene_name` | string | 是 | Gene symbol (e.g. "MYC", "TP53"). |
| `assay_ids` | array of string | 是 | List of assay identifiers. |
| `region` | string \| null | 否 | Input region as "chr1:1000000-1393216". If omitted, auto-centered on the variant position. 默认 `None` |

??? note "上游说明（英文原文）"

    Predict how a variant affects expression of a nearby gene.

    Uses CAGE tracks with TSS-windowed-max and RNA tracks with exon-sum
    quantification, then computes fold change vs reference.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'predict_variant_effect_on_gene',
    {'oracle_name': '...',
     'position': '...',
     'ref_allele': '...',
     'alt_alleles': [],
     'gene_name': '...',
     'assay_ids': []},
)
```

## `score_ism` {#score_ism}

Performs in-silico saturation mutagenesis in a window around a position, returning per-base importance for motif logos.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `center`, `assay_ids`
- **关键词：**ism, saturation mutagenesis, base importance, motif, attribution

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `center` | string | 是 | Variant / motif centre as ``"chrom:pos"`` (1-based). |
| `assay_ids` | array of string | 是 | Track id(s) to score; the first drives the importance profile. |
| `window` | integer | 否 | Window size in bp (default 25). 默认 `25` |
| `genome` | string \| null | 否 | Reference FASTA path; defaults to the oracle's ``reference_fasta``. 默认 `None` |

??? note "上游说明（英文原文）"

    In-silico saturation mutagenesis (ISM) around a variant.

    Sweeps every single-base substitution in a ``window``-bp window centred on
    ``center`` and scores the variant effect on ``assay_ids\[0\]``, returning a
    per-position importance profile (the reference base's disruption) suitable
    for a motif logo. Reveals which bases the oracle actually reads — e.g. a
    disrupted TF motif. Works with any loaded oracle (AlphaGenome, ChromBPNet,
    LegNet, Borzoi, EPInformer-seq, ...).

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'score_ism',
    {'oracle_name': '...', 'center': '...', 'assay_ids': []},
)
```

## `score_prediction_region` {#score_prediction_region}

Predicts a region and summarises signal within a specified sub-region using mean, max, sum or median.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `region`, `assay_ids`, `score_region`
- **关键词：**sub-region scoring, peak, promoter, sequence-to-function

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `region` | string | 是 | Input region as "chr1:1000000-1393216". |
| `assay_ids` | array of string | 是 | List of assay identifiers (e.g. \["DNASE:K562"\]). |
| `score_region` | string | 是 | Sub-region to score as "chr1:1050000-1051000". |
| `scoring_strategy` | string | 否 | How to summarise bins — mean, max, sum, or median. 默认 `mean` |

??? note "上游说明（英文原文）"

    Predict for a region and score a sub-region within the output window.

    Useful for quantifying signal at a specific peak, promoter, or element
    rather than the full output window.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'score_prediction_region',
    {'oracle_name': '...',
     'region': '...',
     'assay_ids': [],
     'score_region': '...'},
)
```

## `score_variant_batch` {#score_variant_batch}

Scores multiple variants through multi-layer analysis and returns a table ranked by effect magnitude.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `variants`, `assay_ids`
- **关键词：**batch variants, vcf, variant ranking, regulatory effect

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `variants` | array of object | 是 | List of variant dicts (schema above). |
| `assay_ids` | array of string | 是 | Track identifiers to score. Pass an empty list to let the oracle score all available tracks (recommended for AlphaGenome / Enformer / Borzoi). |
| `gene_name` | string \| null | 否 | Optional gene for expression scoring. 默认 `None` |
| `top_n` | integer | 否 | Return only the top N variants by effect (default 20). 默认 `20` |
| `user_prompt` | string \| null | 否 | Original user prompt, rendered at the top of the report. 默认 `None` |

??? note "上游说明（英文原文）"

    Score a batch of variants and rank by effect magnitude.

    Processes multiple variants through multi-layer analysis and returns a
    ranked table. Claude can parse VCF content (or any free-text variant
    list) and construct the ``variants`` argument from it.

    **Variant dict schema** — each entry in ``variants`` must be a dict with:

    - ``chrom`` (str): chromosome, e.g. ``"chr1"``
    - ``pos``   (int): 1-based genomic coordinate
    - ``ref``   (str): reference allele, e.g. ``"G"``
    - ``alt``   (str): alternate allele, e.g. ``"T"``
    - ``id``    (str, optional): label, e.g. ``"rs12740374"`` — defaults to
    ``"chrom:pos_ref&gt;alt"`` if omitted

    Example::

    variants = \[
    \{"chrom": "chr1", "pos": 109274968, "ref": "G", "alt": "T", "id": "rs12740374"\},
    \{"chrom": "chr1", "pos": 109275684, "ref": "G", "alt": "T", "id": "rs1626484"\},
    \]

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'score_variant_batch',
    {'oracle_name': '...', 'variants': [], 'assay_ids': []},
)
```

## `score_variant_effect_at_region` {#score_variant_effect_at_region}

Predicts a variant effect and scores ref/alt signal in a chosen sub-region or a window around the variant.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `position`, `ref_allele`, `alt_alleles`, `assay_ids`
- **关键词：**variant effect, sub-region, abs_max, scoring strategy

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `position` | string | 是 | Variant position as "chr1:1050000". |
| `ref_allele` | string | 是 | Reference allele. |
| `alt_alleles` | array of string | 是 | Alternate alleles. |
| `assay_ids` | array of string | 是 | List of assay identifiers. |
| `region` | string \| null | 否 | Input region as "chr1:1000000-1393216". If omitted, auto-centered on the variant position. 默认 `None` |
| `score_region` | string \| null | 否 | Sub-region to score (e.g. "chr1:1050000-1051000"). 默认 `None` |
| `at_variant` | boolean | 否 | If true, score around the variant position instead. 默认 `False` |
| `window_bins` | integer | 否 | Bins on each side when at_variant is true (default 1). 默认 `1` |
| `scoring_strategy` | string | 否 | mean, max, sum, median, or abs_max. 默认 `mean` |

??? note "上游说明（英文原文）"

    Predict a variant effect and score it at a specific region or the variant site.

    Two modes:
    - score_region="chr1:X-Y": score ref/alt in that sub-region.
    - at_variant=True: score ref/alt in a window around the variant position.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'score_variant_effect_at_region',
    {'oracle_name': '...',
     'position': '...',
     'ref_allele': '...',
     'alt_alleles': [],
     'assay_ids': []},
)
```

## `simulate_integration` {#simulate_integration}

Simulates inserting a construct sequence at a genomic position and scores disruption across regulatory layers.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`oracle_name`, `position`, `construct_sequence`, `assay_ids`
- **关键词：**transgene, viral vector, integration site, aavs1, construct

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `oracle_name` | string | 是 | A loaded oracle name. |
| `position` | string | 是 | Insertion point as "chr1:1050000". |
| `construct_sequence` | string | 是 | DNA sequence to insert. |
| `assay_ids` | array of string | 是 | List of assay identifiers for multi-layer scoring. |
| `gene_name` | string \| null | 否 | Optional gene for expression scoring. 默认 `None` |
| `description` | string \| null | 否 | Optional short label (e.g. "AAV-GFP at AAVS1"). 默认 `None` |
| `user_prompt` | string \| null | 否 | Original user prompt, rendered at the top of the report. 默认 `None` |

??? note "上游说明（英文原文）"

    Simulate inserting a construct at a genomic position and score disruption.

    Compares wild-type vs insertion predictions across all regulatory layers.
    Predicts how a viral vector, transgene cassette, or other construct would
    affect local chromatin, TF binding, and gene expression.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'chorus',
    'simulate_integration',
    {'oracle_name': '...',
     'position': '...',
     'construct_sequence': '...',
     'assay_ids': []},
)
```
