# AlphaGenome MCP (JoLab) (`alphagenome-jolab`)

DNAKit 的 MCP 工具索引收录了这个服务器的 **23** 个 DNA/基因相关工具（服务器共列出 24 个，2026-10-16 读取）。工具在用户自己安装服务器后运行，DNAKit 本身不包含它们。调用方式见 [调用索引中的工具](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found)。下面的工具说明是上游英文原文。

| 项目 | 内容 |
| --- | --- |
| 源码 | [taehojo/alphagenome-mcp @ `0c9f2b8`](https://github.com/taehojo/alphagenome-mcp/tree/0c9f2b82c7c2f4cda7e6ba615b4d96dca61645dc) |
| 许可证 | MIT |
| 启动命令 | `alphagenome-mcp` |
| 传输方式 | stdio |
| 安装 | Install with npm install -g @jolab/alphagenome-mcp (Node.js &gt;=18, installed executable alphagenome-mcp) plus pip install 'alphagenome&gt;=0.9.0' numpy in a Python &gt;=3.10 interpreter; set ALPHAGENOME_API_KEY (and optionally ALPHAGENOME_PYTHON to that interpreter, ALPHAGENOME_TIMEOUT_MS). |
| 必需环境变量 | `ALPHAGENOME_API_KEY`（保密） |
| 副作用 | Sends variant coordinates and alleles through a Python subprocess to the Google DeepMind AlphaGenome API (Atlas lookups or live inference) using the caller's API key; no local files are written by the tools. The upstream README states the API is free for non-commercial use only, so commercial use is restricted and requires accepting DeepMind's AlphaGenome Terms of Use and rate limits. |

| 工具 | 类别 | 简介 |
| --- | --- | --- |
| [`analyze_gwas_locus`](#analyze_gwas_locus) | `prediction_model` | Scores 1-100 variants of a locus with live AlphaGenome inference and ranks them by the largest absolute quantile across modalities. |
| [`annotate_regulatory_context`](#annotate_regulatory_context) | `prediction_model` | Runs live AlphaGenome inference for one hg38 variant and reports predicted effects across accessibility, histone, TF binding, CAGE, RNA-seq, splicing and conta… |
| [`assess_pathogenicity`](#assess_pathogenicity) | `prediction_model` | Returns the strongest predicted effect per scorer, the largest quantile and, for SNVs from the Atlas, the AVI score of one hg38 variant; no classification. |
| [`atlas_lookup_variant`](#atlas_lookup_variant) | `prediction_model` | Looks up precomputed AlphaGenome Atlas scores for one hg38 single-nucleotide variant, returning top tracks per scorer with quantile, gene, tissue and assay. |
| [`atlas_lookup_variants`](#atlas_lookup_variants) | `prediction_model` | Looks up precomputed AlphaGenome Atlas scores for up to 500 hg38 SNVs and ranks them, listing variants not found or rejected. |
| [`atlas_scan_region`](#atlas_scan_region) | `prediction_model` | Ranks every possible single-nucleotide substitution in an hg38 region (up to 10 kb, or 50 kb) using precomputed AlphaGenome Atlas scores. |
| [`batch_modality_screen`](#batch_modality_screen) | `prediction_model` | Ranks up to 100 variants by live AlphaGenome-predicted effect within one modality: expression, splicing, TF binding or chromatin. |
| [`batch_pathogenicity_filter`](#batch_pathogenicity_filter) | `prediction_model` | Keeps and ranks variants whose AlphaGenome-predicted effect reaches a calibrated quantile threshold (default 0.99); no pathogenicity classification. |
| [`batch_score_variants`](#batch_score_variants) | `prediction_model` | Scores up to 100 variants with AlphaGenome (SNVs via Atlas, others live) and ranks them by RNA-seq, splice or regulatory-impact metric. |
| [`batch_tissue_comparison`](#batch_tissue_comparison) | `prediction_model` | Ranks up to 100 variants by live AlphaGenome-predicted effect separately within each of up to 10 tissues. |
| [`compare_alleles`](#compare_alleles) | `prediction_model` | Ranks up to 20 alternate alleles at one hg38 position by live AlphaGenome-predicted effect. |
| [`compare_protective_risk`](#compare_protective_risk) | `prediction_model` | Compares AlphaGenome-predicted effect sizes per modality for two variants labeled protective and risk by the caller. |
| [`compare_variants`](#compare_variants) | `prediction_model` | Compares two variants side by side by their strongest AlphaGenome-predicted effect per modality and which is larger per scorer. |
| [`compare_variants_same_gene`](#compare_variants_same_gene) | `prediction_model` | Ranks up to 100 variants by live AlphaGenome-predicted effect, optionally restricting RNA-seq and splice scorers to one gene. |
| [`explain_variant_impact`](#explain_variant_impact) | `prediction_model` | Restates one variant's AlphaGenome-predicted effects (AVI score contributions and strongest effect per modality) as plain sentences. |
| [`generate_variant_report`](#generate_variant_report) | `prediction_model` | Generates a research report of one variant's AlphaGenome-predicted effects with more rows per scorer and AVI feature attributions. |
| [`predict_allele_specific_effects`](#predict_allele_specific_effects) | `prediction_model` | Predicts per-gene, per-tissue expression log fold change of alt vs ref allele (RNA_SEQ) plus expression level with live AlphaGenome inference. |
| [`predict_chromatin_impact`](#predict_chromatin_impact) | `prediction_model` | Predicts ATAC-seq and DNase-seq chromatin accessibility effects of one variant with live AlphaGenome inference. |
| [`predict_expression_impact`](#predict_expression_impact) | `prediction_model` | Predicts RNA-seq log fold change per gene and tissue and CAGE effects of one variant with live AlphaGenome inference. |
| [`predict_splice_impact`](#predict_splice_impact) | `prediction_model` | Predicts splice site, splice site usage and splice junction effects of one variant, with gene and junction, using AlphaGenome. |
| [`predict_tf_binding_impact`](#predict_tf_binding_impact) | `prediction_model` | Predicts transcription factor ChIP-seq binding effects of one variant, with factor and cell type, using AlphaGenome. |
| [`predict_tissue_specific`](#predict_tissue_specific) | `prediction_model` | Reports the strongest AlphaGenome-predicted effect of one variant in each of several requested tissues. |
| [`predict_variant_effect`](#predict_variant_effect) | `prediction_model` | Predicts one variant's regulatory effects per modality (expression, TSS, accessibility, histone, TF binding, splicing) with AlphaGenome, plus AVI from the Atla… |

## `analyze_gwas_locus` {#analyze_gwas_locus}

Scores 1-100 variants of a locus with live AlphaGenome inference and ranks them by the largest absolute quantile across modalities.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`variants`
- **关键词：**gwas locus, fine-mapping, variant prioritization, alphagenome, indel, quantile

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variants` | array of object | 是 | Variants to score (1-100) |
| `chromosome` | string | 否 | Optional: locus chromosome, for the label |
| `start` | number | 否 | Optional: locus start, for the label |
| `end` | number | 否 | Optional: locus end, for the label |

??? note "上游说明（英文原文）"

    Rank the variants of a locus by predicted effect (largest absolute quantile across modalities), to prioritize candidates for follow-up. For single-nucleotide variants only, atlas_lookup_variants or atlas_scan_region is faster and adds the AVI score.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'analyze_gwas_locus',
    {'variants': []},
)
```

## `annotate_regulatory_context` {#annotate_regulatory_context}

Runs live AlphaGenome inference for one hg38 variant and reports predicted effects across accessibility, histone, TF binding, CAGE, RNA-seq, splicing and contact maps.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alt`
- **关键词：**alphagenome, regulatory variant, chromatin accessibility, histone marks, tf binding, splicing, hg38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | 是 | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | 否 | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |

??? note "上游说明（英文原文）"

    The predicted effects of a variant across every regulatory modality at once: accessibility, histone marks, TF binding, CAGE, RNA-seq, splice sites, polyadenylation and contact maps. Shows where the predicted effect concentrates; it does not label the variant.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'annotate_regulatory_context',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```

## `assess_pathogenicity` {#assess_pathogenicity}

Returns the strongest predicted effect per scorer, the largest quantile and, for SNVs from the Atlas, the AVI score of one hg38 variant; no classification.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alt`
- **关键词：**alphagenome, avi score, variant effect, atlas, snv, indel, hg38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | 是 | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | 否 | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |
| `source` | string | 否 | Optional: where the answer comes from (default: auto). auto = the precomputed AlphaGenome Atlas for single-nucleotide substitutions, live inference for everything else (indels, multi-nucleotide variants); falls back to live only when the Atlas does not hold the variant. atlas = Atlas only, errors instead of falling back. live = always run the model. The result always states which source answered. 可选值：`auto`, `atlas`, `live` |
| `scorers` | array of string | 否 | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |

??? note "上游说明（英文原文）"

    Predicted effect size of a variant across modalities, for prioritization. The tool name is kept for compatibility: it does NOT classify a variant as pathogenic or benign, and `classification` is always null.

    Returns the strongest effect per scorer (score, calibrated quantile, where it was seen), the largest absolute quantile, and, for a single-nucleotide variant answered from the Atlas, the AlphaGenome Variant Impact (AVI) score. `avi_score` is null on the live path, because the AVI score is served by the Atlas only.

    Source: a single-nucleotide variant is answered from the precomputed AlphaGenome Atlas; an indel or multi-nucleotide variant runs live inference (score_variant). Both return the same scorers in the same shape. Chosen automatically, overridable with `source`, and always stated in the result.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "How large is the predicted effect of chr19:44908684 T&gt;C?"

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'assess_pathogenicity',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```

## `atlas_lookup_variant` {#atlas_lookup_variant}

Looks up precomputed AlphaGenome Atlas scores for one hg38 single-nucleotide variant, returning top tracks per scorer with quantile, gene, tissue and assay.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alt`
- **关键词：**alphagenome, atlas, snv, avi score, precomputed scores, hg38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference base (one of A, C, G, T). Must match hg38 at this position. |
| `alt` | string | 是 | Alternate base (one of A, C, G, T) |
| `scorers` | array of string | 否 | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |
| `top_n` | number | 否 | Rows to return (default: 25, max: 100) |

??? note "上游说明（英文原文）"

    Look up the precomputed AlphaGenome scores of one single-nucleotide variant. No model call, so it answers in seconds.

    Returns, per scorer, the strongest tracks for the variant ranked by absolute score, each with its calibrated quantile, gene, tissue or cell type, and assay. Default scorers: AVI_SCORE plus one per modality (RNA_SEQ, CAGE, DNASE, CHIP_HISTONE, CHIP_TF, SPLICE_SITES).

    If the reference base does not match hg38, the Atlas says which base it expected and that message is returned as a validation error.

    The Atlas holds precomputed AlphaGenome scores for single-nucleotide substitutions on the human reference genome (hg38, chr1-22, chrX, chrY). Indels and multi-nucleotide variants are not in it; use predict_variant_effect for those.

    The response is a summary, never a full score matrix: ranked rows only, capped at top_n (default 25, max 100) and at 40,000 characters.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Look up chr19:44908684 T&gt;C in the AlphaGenome Atlas"

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'atlas_lookup_variant',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```

## `atlas_lookup_variants` {#atlas_lookup_variants}

Looks up precomputed AlphaGenome Atlas scores for up to 500 hg38 SNVs and ranks them, listing variants not found or rejected.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`variants`
- **关键词：**alphagenome, atlas, batch snv, avi score, variant ranking, hg38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variants` | array of object | 是 | Single-nucleotide variants to look up (1-500) |
| `scorers` | array of string | 否 | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |
| `top_n` | number | 否 | Rows to return (default: 25, max: 100) |

??? note "上游说明（英文原文）"

    Look up precomputed AlphaGenome scores for up to 500 single-nucleotide variants in one call and rank them.

    Returns one row per variant (its strongest score and where it was seen), ranked. Default scorer: AVI_SCORE (AlphaGenome Variant Impact), one number per variant. With several scorers the ranking uses the largest absolute quantile. Variants the Atlas does not hold and variants it rejects (for example a reference base that does not match hg38) are listed separately with the reason; they do not fail the call.

    The Atlas holds precomputed AlphaGenome scores for single-nucleotide substitutions on the human reference genome (hg38, chr1-22, chrX, chrY). Indels and multi-nucleotide variants are not in it; use predict_variant_effect for those.

    The response is a summary, never a full score matrix: ranked rows only, capped at top_n (default 25, max 100) and at 40,000 characters.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Rank these 200 GWAS SNPs by their Atlas scores"

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'atlas_lookup_variants',
    {'variants': []},
)
```

## `atlas_scan_region` {#atlas_scan_region}

Ranks every possible single-nucleotide substitution in an hg38 region (up to 10 kb, or 50 kb) using precomputed AlphaGenome Atlas scores.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `start`, `end`
- **关键词：**alphagenome, saturation mutagenesis, region scan, atlas, avi score, hg38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `start` | number | 是 | Start position (1-based, hg38) |
| `end` | number | 是 | End position (greater than start; at most start + 10,000, or start + 50,000 with allow_large_region) |
| `allow_large_region` | boolean | 否 | Set to true to scan more than 10,000 bp (up to 50,000 bp). Slower, and the result may be incomplete (default: false) |
| `scorers` | array of string | 否 | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |
| `top_n` | number | 否 | Rows to return (default: 25, max: 100) |

??? note "上游说明（英文原文）"

    Scan a genomic region in the AlphaGenome Atlas: every possible single-nucleotide substitution in the interval, ranked. Answers "which positions in this region matter most?" without running the model.

    Region width: at most 10,000 bp. Up to 50,000 bp only with allow_large_region=true. A scan is one API request per 32 bp under a requests-per-minute quota, so a large scan can take minutes. If the quota or the time limit stops a scan early, the partial result is returned, marked "Incomplete", with the range that was really scanned; the ranking then covers that range only.

    Default scorer: AVI_SCORE (about 7 seconds for 2,000 bp). Multi-track scorers are much slower (2,000 bp: DNASE 17 s, CHIP_TF 86 s). Scorers with one row per gene or junction (RNA_SEQ, SPLICE_JUNCTIONS, ...) cannot be used for a scan: scan with AVI_SCORE, then use atlas_lookup_variant on the top variants.

    The Atlas holds precomputed AlphaGenome scores for single-nucleotide substitutions on the human reference genome (hg38, chr1-22, chrX, chrY). Indels and multi-nucleotide variants are not in it; use predict_variant_effect for those.

    The response is a summary, never a full score matrix: ranked rows only, capped at top_n (default 25, max 100) and at 40,000 characters.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Scan chr17:49209289-49211289 and show the 10 substitutions with the largest predicted effect"

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'atlas_scan_region',
    {'chromosome': '...', 'start': 0.0, 'end': 0.0},
)
```

## `batch_modality_screen` {#batch_modality_screen}

Ranks up to 100 variants by live AlphaGenome-predicted effect within one modality: expression, splicing, TF binding or chromatin.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`variants`, `modality`
- **关键词：**alphagenome, variant ranking, expression, splicing, tf binding, chromatin

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variants` | array of object | 是 | Variants to score (1-100) |
| `modality` | string | 是 | 可选值：`expression`, `splicing`, `tf_binding`, `chromatin` |

??? note "上游说明（英文原文）"

    Rank several variants by predicted effect within one modality: expression (RNA_SEQ, CAGE), splicing (SPLICE_SITES, SPLICE_SITE_USAGE), tf_binding (CHIP_TF) or chromatin (DNASE, ATAC).

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'batch_modality_screen',
    {'variants': [], 'modality': 'expression'},
)
```

## `batch_pathogenicity_filter` {#batch_pathogenicity_filter}

Keeps and ranks variants whose AlphaGenome-predicted effect reaches a calibrated quantile threshold (default 0.99); no pathogenicity classification.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`variants`
- **关键词：**alphagenome, variant filter, quantile threshold, avi score, prioritization

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variants` | array of object | 是 | Variants to score (1-100) |
| `threshold` | number | 否 | Smallest absolute quantile to keep, between 0 and 1 (default: 0.99) |
| `source` | string | 否 | Optional: where the answer comes from (default: auto). auto = the precomputed AlphaGenome Atlas for single-nucleotide substitutions, live inference for everything else (indels, multi-nucleotide variants); falls back to live only when the Atlas does not hold the variant. atlas = Atlas only, errors instead of falling back. live = always run the model. The result always states which source answered. 可选值：`auto`, `atlas`, `live` |
| `scorers` | array of string | 否 | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |

??? note "上游说明（英文原文）"

    Keep the variants whose predicted effect reaches a threshold, ranked. The tool name is kept for compatibility: it filters on predicted effect size, NOT on pathogenicity, and classifies nothing.

    `threshold` is the smallest absolute calibrated quantile (0 to 1) a variant must reach to be kept (default 0.99): the AVI score's quantile for variants answered from the Atlas, the largest quantile across modalities for live inference. Variants are routed per variant and reported per source; groups from different sources are not comparable.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Which of these variants have a predicted effect above the 99.9th percentile?"

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'batch_pathogenicity_filter',
    {'variants': []},
)
```

## `batch_score_variants` {#batch_score_variants}

Scores up to 100 variants with AlphaGenome (SNVs via Atlas, others live) and ranks them by RNA-seq, splice or regulatory-impact metric.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`variants`, `scoring_metric`
- **关键词：**alphagenome, batch scoring, variant ranking, atlas, rna-seq, splicing

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variants` | array of object | 是 | Variants to score (1-100) |
| `scoring_metric` | string | 是 | What to rank by 可选值：`rna_seq`, `splice`, `regulatory_impact`, `combined` |
| `top_n` | number | 否 | Variants to return per group (default: 10, max: 100) |
| `source` | string | 否 | Optional: where the answer comes from (default: auto). auto = the precomputed AlphaGenome Atlas for single-nucleotide substitutions, live inference for everything else (indels, multi-nucleotide variants); falls back to live only when the Atlas does not hold the variant. atlas = Atlas only, errors instead of falling back. live = always run the model. The result always states which source answered. 可选值：`auto`, `atlas`, `live` |
| `scorers` | array of string | 否 | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |

??? note "上游说明（英文原文）"

    Score up to 100 variants and rank them by predicted effect.

    Each variant is routed on its own: single-nucleotide variants to the Atlas, the rest to live inference. A mixed batch comes back as two separately ranked groups, because the Atlas group is ranked by the AVI score and live inference has no AVI score; the two must not be compared. The result reports how many variants came from each source and how many fell back.

    Scoring metric (used when `scorers` is not given): rna_seq = RNA_SEQ, splice = SPLICE_SITES, regulatory_impact and combined = AVI_SCORE from the Atlas and every modality (ranked by the largest absolute quantile) from live inference.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Score these 50 variants and show me the top 10 by predicted effect"

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'batch_score_variants',
    {'variants': [], 'scoring_metric': 'rna_seq'},
)
```

## `batch_tissue_comparison` {#batch_tissue_comparison}

Ranks up to 100 variants by live AlphaGenome-predicted effect separately within each of up to 10 tissues.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`variants`, `tissues`
- **关键词：**alphagenome, tissue-specific, variant ranking, uberon, tissues

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variants` | array of object | 是 | Variants to score (1-100) |
| `tissues` | array of string | 是 | Tissue names or ontology CURIEs (1-10) |

??? note "上游说明（英文原文）"

    Rank several variants by predicted effect within each of several tissues: one ranking per tissue.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'batch_tissue_comparison',
    {'variants': [], 'tissues': []},
)
```

## `compare_alleles` {#compare_alleles}

Ranks up to 20 alternate alleles at one hg38 position by live AlphaGenome-predicted effect.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alts`
- **关键词：**alphagenome, alternate alleles, allele comparison, hg38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference allele (A, C, G, T; more than one base for an indel) |
| `alts` | array of string | 是 | Alternate alleles (1-20) |

??? note "上游说明（英文原文）"

    Rank the alternate alleles of one position by predicted effect.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'compare_alleles',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alts': []},
)
```

## `compare_protective_risk` {#compare_protective_risk}

Compares AlphaGenome-predicted effect sizes per modality for two variants labeled protective and risk by the caller.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`protective_variant`, `risk_variant`
- **关键词：**alphagenome, protective allele, risk allele, variant comparison

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `protective_variant` | object | 是 |  |
| `risk_variant` | object | 是 |  |

??? note "上游说明（英文原文）"

    Two variants side by side, labelled as the caller names them. The tool compares predicted effect sizes per modality; it does not judge which allele is protective or a risk.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'compare_protective_risk',
    {'protective_variant': {}, 'risk_variant': {}},
)
```

## `compare_variants` {#compare_variants}

Compares two variants side by side by their strongest AlphaGenome-predicted effect per modality and which is larger per scorer.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`variant1`, `variant2`
- **关键词：**alphagenome, variant comparison, rsid, modality, quantile

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variant1` | object | 是 |  |
| `variant2` | object | 是 |  |

??? note "上游说明（英文原文）"

    Two variants side by side: the strongest predicted effect of each modality for both, and which of the two has the larger absolute quantile per scorer. A comparison of predicted effect sizes, not of severity.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Compare APOE rs429358 and rs7412"

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'compare_variants',
    {'variant1': {}, 'variant2': {}},
)
```

## `compare_variants_same_gene` {#compare_variants_same_gene}

Ranks up to 100 variants by live AlphaGenome-predicted effect, optionally restricting RNA-seq and splice scorers to one gene.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`variants`
- **关键词：**alphagenome, gene symbol, variant ranking, splicing, expression

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `variants` | array of object | 是 | Variants to score (1-100) |
| `gene_name` | string | 否 | Optional: gene symbol (e.g., APOE) |

??? note "上游说明（英文原文）"

    Rank several variants by their predicted effect on one gene. With gene_name, the gene-level scorers (RNA_SEQ, SPLICE_SITES) are restricted to that gene.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'compare_variants_same_gene',
    {'variants': []},
)
```

## `explain_variant_impact` {#explain_variant_impact}

Restates one variant's AlphaGenome-predicted effects (AVI score contributions and strongest effect per modality) as plain sentences.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alt`
- **关键词：**alphagenome, variant interpretation, avi score, natural language summary, hg38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | 是 | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | 否 | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |
| `source` | string | 否 | Optional: where the answer comes from (default: auto). auto = the precomputed AlphaGenome Atlas for single-nucleotide substitutions, live inference for everything else (indels, multi-nucleotide variants); falls back to live only when the Atlas does not hold the variant. atlas = Atlas only, errors instead of falling back. live = always run the model. The result always states which source answered. 可选值：`auto`, `atlas`, `live` |
| `scorers` | array of string | 否 | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |

??? note "上游说明（英文原文）"

    Plain sentences that restate a variant's predicted effects: the AVI score and its largest contributions (from the Atlas), then the strongest effect of each modality ordered by absolute quantile, with the direction for signed scorers. Descriptive only: the sentences restate returned numbers and make no statement about pathogenicity.

    Source: a single-nucleotide variant is answered from the precomputed AlphaGenome Atlas; an indel or multi-nucleotide variant runs live inference (score_variant). Both return the same scorers in the same shape. Chosen automatically, overridable with `source`, and always stated in the result.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Explain the predicted effect of chr17:49210289 C&gt;T in plain language"

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'explain_variant_impact',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```

## `generate_variant_report` {#generate_variant_report}

Generates a research report of one variant's AlphaGenome-predicted effects with more rows per scorer and AVI feature attributions.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alt`
- **关键词：**alphagenome, variant report, avi score, feature importance, hg38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | 是 | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | 否 | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |
| `source` | string | 否 | Optional: where the answer comes from (default: auto). auto = the precomputed AlphaGenome Atlas for single-nucleotide substitutions, live inference for everything else (indels, multi-nucleotide variants); falls back to live only when the Atlas does not hold the variant. atlas = Atlas only, errors instead of falling back. live = always run the model. The result always states which source answered. 可选值：`auto`, `atlas`, `live` |
| `scorers` | array of string | 否 | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |

??? note "上游说明（英文原文）"

    A fuller report of one variant's predicted molecular effects: more rows per scorer than predict_variant_effect and, from the Atlas, the AVI score with its feature attributions (AVI_SCORE_FEATURE_IMPORTANCE). It is a research summary, not a clinical report: it contains no pathogenicity classification and no recommendation.

    Source: a single-nucleotide variant is answered from the precomputed AlphaGenome Atlas; an indel or multi-nucleotide variant runs live inference (score_variant). Both return the same scorers in the same shape. Chosen automatically, overridable with `source`, and always stated in the result.

    The response is a summary, never a full score matrix: ranked rows only, capped at top_n (default 25, max 100) and at 40,000 characters.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Generate a report for chr19:44908684 T&gt;C"

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'generate_variant_report',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```

## `predict_allele_specific_effects` {#predict_allele_specific_effects}

Predicts per-gene, per-tissue expression log fold change of alt vs ref allele (RNA_SEQ) plus expression level with live AlphaGenome inference.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alt`
- **关键词：**alphagenome, allele-specific expression, log fold change, rna-seq, eqtl-like

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | 是 | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | 否 | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |

??? note "上游说明（英文原文）"

    Predicted expression with the alternate allele against the reference allele: the RNA_SEQ scorer is that log fold change per gene and tissue, and RNA_SEQ_ACTIVE gives the expression level alongside it.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'predict_allele_specific_effects',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```

## `predict_chromatin_impact` {#predict_chromatin_impact}

Predicts ATAC-seq and DNase-seq chromatin accessibility effects of one variant with live AlphaGenome inference.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alt`
- **关键词：**alphagenome, atac-seq, dnase-seq, chromatin accessibility, regulatory variant

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | 是 | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | 否 | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |

??? note "上游说明（英文原文）"

    Predicted chromatin accessibility effects of a variant (ATAC-seq and DNase-seq).

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'predict_chromatin_impact',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```

## `predict_expression_impact` {#predict_expression_impact}

Predicts RNA-seq log fold change per gene and tissue and CAGE effects of one variant with live AlphaGenome inference.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alt`
- **关键词：**alphagenome, gene expression, rna-seq, cage, log fold change

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | 是 | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | 否 | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |

??? note "上游说明（英文原文）"

    Predicted gene expression effects of a variant: RNA-seq (log fold change per gene and tissue) and CAGE.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'predict_expression_impact',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```

## `predict_splice_impact` {#predict_splice_impact}

Predicts splice site, splice site usage and splice junction effects of one variant, with gene and junction, using AlphaGenome.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alt`
- **关键词：**alphagenome, splicing, splice junction, splice site usage, sqtl

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | 是 | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | 否 | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |

??? note "上游说明（英文原文）"

    Predicted splicing effects of a variant: splice sites, splice site usage and splice junctions, with the gene and junction of each.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'predict_splice_impact',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```

## `predict_tf_binding_impact` {#predict_tf_binding_impact}

Predicts transcription factor ChIP-seq binding effects of one variant, with factor and cell type, using AlphaGenome.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alt`
- **关键词：**alphagenome, tf chip-seq, transcription factor, binding site, regulatory variant

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | 是 | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | 否 | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |

??? note "上游说明（英文原文）"

    Predicted transcription factor binding effects of a variant (TF ChIP-seq), with the factor and cell type of each.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'predict_tf_binding_impact',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```

## `predict_tissue_specific` {#predict_tissue_specific}

Reports the strongest AlphaGenome-predicted effect of one variant in each of several requested tissues.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alt`
- **关键词：**alphagenome, tissue-specific, brain, liver, heart, variant effect

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | 是 | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissues` | array of string | 否 | Tissue names (brain, neuron, blood, liver, heart, lung, kidney) or ontology CURIEs. Default: brain, liver, heart |

??? note "上游说明（英文原文）"

    The strongest predicted effect of a variant in each of several tissues: the same scores, filtered to the tracks of one tissue at a time.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Compare the predicted effect of chr19:44908684 T&gt;C in brain, liver and heart"

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'predict_tissue_specific',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```

## `predict_variant_effect` {#predict_variant_effect}

Predicts one variant's regulatory effects per modality (expression, TSS, accessibility, histone, TF binding, splicing) with AlphaGenome, plus AVI from the Atlas.

- **类别：**`prediction_model`（DNA sequence models and variant-effect prediction models）
- **必填参数：**`chromosome`, `position`, `ref`, `alt`
- **关键词：**alphagenome, variant effect prediction, avi score, regulatory, splicing, hg38

| 参数 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `chromosome` | string | 是 | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | 是 | Genomic position (1-based, hg38) |
| `ref` | string | 是 | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | 是 | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | 否 | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |
| `output_types` | array of string | 否 | Optional: modalities to report (default: one scorer per modality) |
| `source` | string | 否 | Optional: where the answer comes from (default: auto). auto = the precomputed AlphaGenome Atlas for single-nucleotide substitutions, live inference for everything else (indels, multi-nucleotide variants); falls back to live only when the Atlas does not hold the variant. atlas = Atlas only, errors instead of falling back. live = always run the model. The result always states which source answered. 可选值：`auto`, `atlas`, `live` |
| `scorers` | array of string | 否 | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |

??? note "上游说明（英文原文）"

    Predicted regulatory effect of a genetic variant, per modality: the strongest tracks of each scorer (expression, transcription start, chromatin accessibility, histone marks, transcription factor binding, splicing), each with its score, calibrated quantile, gene, tissue or cell type, and assay. From the Atlas the AlphaGenome Variant Impact (AVI) score is included.

    Source: a single-nucleotide variant is answered from the precomputed AlphaGenome Atlas; an indel or multi-nucleotide variant runs live inference (score_variant). Both return the same scorers in the same shape. Chosen automatically, overridable with `source`, and always stated in the result.

    The response is a summary, never a full score matrix: ranked rows only, capped at top_n (default 25, max 100) and at 40,000 characters.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Analyze chr19:44908684 T&gt;C with AlphaGenome"

安装服务器后通过 DNAKit 调用：

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'predict_variant_effect',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```
