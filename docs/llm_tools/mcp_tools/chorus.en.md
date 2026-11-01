# Chorus (`chorus`)

DNAKit's MCP tool index keeps **17** DNA/gene tools of this server (24 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [pinellolab/chorus @ `cc8aabe`](https://github.com/pinellolab/chorus/tree/cc8aabe9f5e4e562f0fc8b8bfa838fb100cba52b) |
| License | MIT |
| Command | `chorus-mcp` |
| Transports | stdio |
| Installation | Install the pinellolab/chorus source project and provision its oracle environments, model weights and references separately. The PyPI package named chorus is unrelated. Set executable_path to the installed chorus-mcp. |
| Effects | Predictions can consume GPU/CPU, write reports, and download missing model assets inside the upstream server. Pre-provision assets when offline. |

| Tool | Category | Summary |
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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `region`, `replacement_sequence`, `assay_ids`
- **Keywords:** promoter swap, enhancer replacement, sequence-to-function, alphagenome, enformer

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `region` | string | yes | Region to replace as "chr1:1000000-1001000". |
| `replacement_sequence` | string | yes | DNA sequence to insert in place of the region. |
| `assay_ids` | array of string | yes | List of assay identifiers for multi-layer scoring. |
| `gene_name` | string \| null | no | Optional gene for expression scoring. Default `None` |
| `description` | string \| null | no | Optional short label for the swap (e.g. "SV40 promoter"). Default `None` |
| `user_prompt` | string \| null | no | Original user prompt, rendered at the top of the report. Default `None` |

??? note "Upstream description"

    Replace a genomic region with a custom sequence and score effects across all layers.

    Compares wild-type vs replacement predictions using the same multi-layer
    scoring as variant analysis (chromatin, TF binding, histone, CAGE, RNA).

    Use cases: promoter swaps, enhancer replacements, regulatory element engineering.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `position`, `ref_allele`, `alt_alleles`, `assay_ids`
- **Keywords:** variant effect prediction, chromatin, tf binding, cage, log2 fold change

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `position` | string | yes | Variant position as ``"chr1:1050000"`` OR an rsID like ``"rs12740374"`` (requires LDlink token; resolves to coords via the LDproxy API). |
| `ref_allele` | string | yes | Reference allele. |
| `alt_alleles` | array of string | yes | Alternate alleles. |
| `assay_ids` | array of string | yes | List of assay identifiers covering different layers (e.g. DNASE, CAGE, ChIP tracks for multi-layer coverage). Pass an empty list or None to score all tracks on oracles that support it (AlphaGenome, Enformer, Borzoi). |
| `gene_name` | string \| null | no | Gene symbol for RNA expression scoring (e.g. "SORT1"). If omitted, the nearest gene is auto-detected. Default `None` |
| `region` | string \| null | no | Input region as "chr1:1000000-1393216". If omitted, auto-centered. Default `None` |
| `igv_raw` | boolean | no | When True, the IGV browser in the HTML report shows raw signal with autoscale instead of the layer-aware rescaled view. Table scores are unaffected. Default `False` |
| `ldlink_token` | string \| null | no | LDlink API token (only used when ``position`` is an rsID). Register free at https://ldlink.nih.gov/?tab=apiaccess or set ``LDLINK_TOKEN``. Default `None` |
| `genome_build` | string | no | Reference build for the rsID lookup (``"grch38"`` / ``"hg38"`` default; or ``"grch37"`` / ``"hg19"``). Only used when ``position`` is an rsID. Default `grch38` |
| `user_prompt` | string \| null | no | The user's original natural-language question. Claude should forward this verbatim whenever calling from an MCP conversation — it is rendered at the top of the report for traceability. Default `None` |

??? note "Upstream description"

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

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `position`, `ref_allele`, `alt_alleles`
- **Keywords:** variant interpretation, track ranking, alphagenome, borzoi, cell types

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `position` | string | yes | Variant position as ``"chr1:109274968"`` OR an rsID like ``"rs12740374"`` (requires LDlink token; resolves to coords via the LDproxy API). |
| `ref_allele` | string | yes | Reference allele (e.g. "G"). |
| `alt_alleles` | array of string | yes | Alternate alleles (e.g. \["T"\]). |
| `gene_name` | string \| null | no | Optional gene for expression analysis. Default `None` |
| `top_n` | integer | no | Number of top tracks per regulatory layer to show. Default `3` |
| `igv_raw` | boolean | no | Default `False` |
| `user_prompt` | string \| null | no | Original user prompt, forwarded into the report header. Default `None` |
| `ranking_metric` | string | no | How to rank tracks and cell types. Default ``"alt_x_abs_effect"`` (``alt_value × \|log2FC\|``) avoids inflating closed-baseline tracks where a small absolute change yields a large fold-change. Use ``"abs_effect"`` for the historical raw \|log2FC\| ranking. Default `alt_x_abs_effect` |
| `min_ref_value` | number | no | Threshold used by ``"abs_effect_min_ref"``. Default `0.0` |
| `ldlink_token` | string \| null | no | LDlink API token (only used when ``position`` is an rsID). Register free at https://ldlink.nih.gov/?tab=apiaccess. Default `None` |
| `genome_build` | string | no | Reference build for the rsID lookup (``"grch38"`` / ``"hg38"`` default; or ``"grch37"`` / ``"hg19"``). Default `grch38` |

??? note "Upstream description"

    Discover which cell types and regulatory layers are most affected by a variant.

    Predicts variant effect across ALL available tracks (thousands for
    Enformer/Borzoi/AlphaGenome, or iterates all models for ChromBPNet/LegNet),
    ranks by ``ranking_metric``, and returns the top hits with a full
    report.

    This is the primary tool for variant interpretation — it tells you WHERE
    the variant has impact without requiring you to pre-select tracks.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `position`, `ref_allele`, `alt_alleles`
- **Keywords:** cell-type discovery, dnase, atac, variant effect, alphagenome

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name (ideally AlphaGenome for broadest cell-type coverage). |
| `position` | string | yes | Variant position as ``"chr1:1050000"`` OR an rsID like ``"rs12740374"`` (requires LDlink token; resolves to coords via the LDproxy API). |
| `ref_allele` | string | yes | Reference allele. |
| `alt_alleles` | array of string | yes | Alternate alleles. |
| `gene_name` | string \| null | no | Optional gene to focus expression analysis on. Default `None` |
| `top_n` | integer | no | Number of top cell types to analyze in detail (default 5). Default `5` |
| `min_effect` | number | no | Minimum \|log2FC\| in DNASE/ATAC to consider a cell type hit (default 0.15). Default `0.15` |
| `user_prompt` | string \| null | no | Original user prompt, forwarded into each sub-report. Default `None` |
| `ranking_metric` | string | no | How to rank cell types. Default ``"alt_x_abs_effect"`` (recommended) ranks by ``alt_signal × \|log2FC\|``, which surfaces the cell type where the post-variant element is most active. Use ``"abs_effect"`` to reproduce the historical \|log2FC\|-only ranking (which has a bias toward closed-baseline cell types — see :func:`chorus.analysis.discovery.discover_cell_types`). Use ``"abs_effect_min_ref"`` with… Default `alt_x_abs_effect` |
| `min_ref_value` | number | no | Threshold for ``"abs_effect_min_ref"`` (ignored otherwise). Default 0 = no filter. Default `0.0` |
| `ldlink_token` | string \| null | no | LDlink API token (only used when ``position`` is an rsID). Default `None` |
| `genome_build` | string | no | Reference build for the rsID lookup (``"grch38"`` / ``"hg38"`` default; or ``"grch37"`` / ``"hg19"``). Default `grch38` |

??? note "Upstream description"

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

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `lead_variant`
- **Keywords:** fine-mapping, gwas, ld proxies, ldlink, causal variant

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `lead_variant` | string | yes | Sentinel variant as "rs12740374" or "chr1:109274968 G&gt;T". |
| `ld_variants` | array \| null | no | Optional list of LD variant dicts (schema above). If omitted, auto-fetched from LDlink. Default `None` |
| `assay_ids` | array \| null | no | Track identifiers. Pass an empty list / None to score all tracks (recommended for AlphaGenome). Default `None` |
| `gene_name` | string \| null | no | Target gene for expression scoring. Default `None` |
| `population` | string | no | 1000 Genomes population for LD lookup (default CEU). Default `CEU` |
| `r2_threshold` | number | no | Minimum r² for LD variants (default 0.8). Default `0.8` |
| `ldlink_token` | string \| null | no | LDlink API token. Register free at https://ldlink.nih.gov/?tab=apiaccess Default `None` |
| `ldlink_timeout` | number | no | HTTP timeout in seconds for the LDlink request (default 30). Increase for slow networks or large LD blocks. Default `30.0` |
| `genome_build` | string | no | Reference build for the LDlink LD lookup. Accepts ``"grch38"`` / ``"hg38"`` (default) or ``"grch37"`` / ``"hg19"``. Default `grch38` |
| `snvs_only` | boolean | no | When True, restrict scoring to single-nucleotide substitutions (drops insertions, deletions, MNVs, and complex multi-base proxies). Default False — indels and multi-allelic LDlink rows are scored by default as of v0.5.5. Set True to reproduce pre-v0.5.5 behavior. Default `False` |
| `user_prompt` | string \| null | no | Original user prompt, rendered at the top of the report. Default `None` |
| `report_top_n` | integer | no | Number of top-ranked variants to render full IGV signal tracks for in the HTML report (default 3). Every variant is scored with the same per-track logic; this only controls how many of the top hits get the (more expensive) interactive browser tracks. Set to 0 to skip the IGV signal tracks entirely. Default `3` |

??? note "Upstream description"

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

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `gene_name`
- **Keywords:** tss, transcription start site, gene annotation, promoter

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `gene_name` | string | yes | Gene symbol (e.g. "GATA1", "TP53"). |

??? note "Upstream description"

    Get transcription start site (TSS) positions for a gene.

Call it through DNAKit once the server is installed:

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

- **Category:** `genome_annotation` (Genes, transcripts, gene IDs, genomic coordinates, assemblies and gene annotation)
- **Required arguments:** `chrom`, `start`, `end`
- **Keywords:** genomic region, gene annotation, chromosome, interval

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chrom` | string | yes | Chromosome (e.g. "chr1"). |
| `start` | integer | yes | Region start position. |
| `end` | integer | yes | Region end position. |

??? note "Upstream description"

    Get gene annotations in a genomic region.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `region`, `assay_ids`
- **Keywords:** sequence-to-function, bedgraph, dnase, cage, genomic region

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `region` | string | yes | Genomic region as "chr1:1000000-1393216". |
| `assay_ids` | array of string | yes | List of assay identifiers (e.g. \["DNASE:K562", "CAGE:K562"\]). |

??? note "Upstream description"

    Make a wild-type prediction for a genomic region.

    Returns per-track summary stats (mean/max/min/std). For small predictions
    the full values are included inline; for large ones a downsampled preview
    is returned and full data is saved as bedgraph files.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `position`, `sequence`, `assay_ids`
- **Keywords:** sequence insertion, in silico, sequence-to-function, regulatory activity

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `position` | string | yes | Insertion point as "chr1:1050000". |
| `sequence` | string | yes | DNA sequence to insert. |
| `assay_ids` | array of string | yes | List of assay identifiers. |

??? note "Upstream description"

    Insert a sequence at a genomic position and predict activity.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `region`, `replacement_sequence`, `assay_ids`
- **Keywords:** region replacement, in silico, sequence-to-function, regulatory activity

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `region` | string | yes | Genomic region to replace as "chr1:1000000-1001000". |
| `replacement_sequence` | string | yes | DNA sequence to insert in place of the region. |
| `assay_ids` | array of string | yes | List of assay identifiers. |

??? note "Upstream description"

    Replace a genomic region with a custom sequence and predict activity.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `position`, `ref_allele`, `alt_alleles`, `assay_ids`
- **Keywords:** variant effect prediction, ref alt, log2 fold change, tracks

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `position` | string | yes | Variant position as "chr1:1050000". |
| `ref_allele` | string | yes | Reference allele (e.g. "A"). |
| `alt_alleles` | array of string | yes | Alternate alleles (e.g. \["G", "T"\]). |
| `assay_ids` | array of string | yes | List of assay identifiers. |
| `region` | string \| null | no | Genomic region as "chr1:1000000-1393216". If omitted, auto-centered on the variant position using the oracle's input window. Default `None` |

??? note "Upstream description"

    Predict the effect of a genetic variant.

    Compares predictions for reference vs alternate alleles and returns
    per-allele effect sizes with summary statistics.

    One forward pass scores the variant across every track in ``assay_ids``
    (or ALL of the oracle's tracks when ``assay_ids`` is omitted), so request
    the full set in a single call — never loop this tool once per track or per
    cell type. To score many variants, use ``score_variant_batch`` or
    ``fine_map_causal_variant`` (each scores a variant exactly once).

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `position`, `ref_allele`, `alt_alleles`, `gene_name`, `assay_ids`
- **Keywords:** gene expression, cage, rna, variant effect, fold change

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `position` | string | yes | Variant position as "chr1:1050000". |
| `ref_allele` | string | yes | Reference allele. |
| `alt_alleles` | array of string | yes | Alternate alleles. |
| `gene_name` | string | yes | Gene symbol (e.g. "MYC", "TP53"). |
| `assay_ids` | array of string | yes | List of assay identifiers. |
| `region` | string \| null | no | Input region as "chr1:1000000-1393216". If omitted, auto-centered on the variant position. Default `None` |

??? note "Upstream description"

    Predict how a variant affects expression of a nearby gene.

    Uses CAGE tracks with TSS-windowed-max and RNA tracks with exon-sum
    quantification, then computes fold change vs reference.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `center`, `assay_ids`
- **Keywords:** ism, saturation mutagenesis, base importance, motif, attribution

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `center` | string | yes | Variant / motif centre as ``"chrom:pos"`` (1-based). |
| `assay_ids` | array of string | yes | Track id(s) to score; the first drives the importance profile. |
| `window` | integer | no | Window size in bp (default 25). Default `25` |
| `genome` | string \| null | no | Reference FASTA path; defaults to the oracle's ``reference_fasta``. Default `None` |

??? note "Upstream description"

    In-silico saturation mutagenesis (ISM) around a variant.

    Sweeps every single-base substitution in a ``window``-bp window centred on
    ``center`` and scores the variant effect on ``assay_ids\[0\]``, returning a
    per-position importance profile (the reference base's disruption) suitable
    for a motif logo. Reveals which bases the oracle actually reads — e.g. a
    disrupted TF motif. Works with any loaded oracle (AlphaGenome, ChromBPNet,
    LegNet, Borzoi, EPInformer-seq, ...).

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `region`, `assay_ids`, `score_region`
- **Keywords:** sub-region scoring, peak, promoter, sequence-to-function

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `region` | string | yes | Input region as "chr1:1000000-1393216". |
| `assay_ids` | array of string | yes | List of assay identifiers (e.g. \["DNASE:K562"\]). |
| `score_region` | string | yes | Sub-region to score as "chr1:1050000-1051000". |
| `scoring_strategy` | string | no | How to summarise bins — mean, max, sum, or median. Default `mean` |

??? note "Upstream description"

    Predict for a region and score a sub-region within the output window.

    Useful for quantifying signal at a specific peak, promoter, or element
    rather than the full output window.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `variants`, `assay_ids`
- **Keywords:** batch variants, vcf, variant ranking, regulatory effect

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `variants` | array of object | yes | List of variant dicts (schema above). |
| `assay_ids` | array of string | yes | Track identifiers to score. Pass an empty list to let the oracle score all available tracks (recommended for AlphaGenome / Enformer / Borzoi). |
| `gene_name` | string \| null | no | Optional gene for expression scoring. Default `None` |
| `top_n` | integer | no | Return only the top N variants by effect (default 20). Default `20` |
| `user_prompt` | string \| null | no | Original user prompt, rendered at the top of the report. Default `None` |

??? note "Upstream description"

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

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `position`, `ref_allele`, `alt_alleles`, `assay_ids`
- **Keywords:** variant effect, sub-region, abs_max, scoring strategy

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `position` | string | yes | Variant position as "chr1:1050000". |
| `ref_allele` | string | yes | Reference allele. |
| `alt_alleles` | array of string | yes | Alternate alleles. |
| `assay_ids` | array of string | yes | List of assay identifiers. |
| `region` | string \| null | no | Input region as "chr1:1000000-1393216". If omitted, auto-centered on the variant position. Default `None` |
| `score_region` | string \| null | no | Sub-region to score (e.g. "chr1:1050000-1051000"). Default `None` |
| `at_variant` | boolean | no | If true, score around the variant position instead. Default `False` |
| `window_bins` | integer | no | Bins on each side when at_variant is true (default 1). Default `1` |
| `scoring_strategy` | string | no | mean, max, sum, median, or abs_max. Default `mean` |

??? note "Upstream description"

    Predict a variant effect and score it at a specific region or the variant site.

    Two modes:
    - score_region="chr1:X-Y": score ref/alt in that sub-region.
    - at_variant=True: score ref/alt in a window around the variant position.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `oracle_name`, `position`, `construct_sequence`, `assay_ids`
- **Keywords:** transgene, viral vector, integration site, aavs1, construct

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `oracle_name` | string | yes | A loaded oracle name. |
| `position` | string | yes | Insertion point as "chr1:1050000". |
| `construct_sequence` | string | yes | DNA sequence to insert. |
| `assay_ids` | array of string | yes | List of assay identifiers for multi-layer scoring. |
| `gene_name` | string \| null | no | Optional gene for expression scoring. Default `None` |
| `description` | string \| null | no | Optional short label (e.g. "AAV-GFP at AAVS1"). Default `None` |
| `user_prompt` | string \| null | no | Original user prompt, rendered at the top of the report. Default `None` |

??? note "Upstream description"

    Simulate inserting a construct at a genomic position and score disruption.

    Compares wild-type vs insertion predictions across all regulatory layers.
    Predicts how a viral vector, transgene cassette, or other construct would
    affect local chromatin, TF binding, and gene expression.

Call it through DNAKit once the server is installed:

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
