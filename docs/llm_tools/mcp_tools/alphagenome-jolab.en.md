# AlphaGenome MCP (JoLab) (`alphagenome-jolab`)

DNAKit's MCP tool index keeps **23** DNA/gene tools of this server (24 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [taehojo/alphagenome-mcp @ `0c9f2b8`](https://github.com/taehojo/alphagenome-mcp/tree/0c9f2b82c7c2f4cda7e6ba615b4d96dca61645dc) |
| License | MIT |
| Command | `alphagenome-mcp` |
| Transports | stdio |
| Installation | Install with npm install -g @jolab/alphagenome-mcp (Node.js &gt;=18, installed executable alphagenome-mcp) plus pip install 'alphagenome&gt;=0.9.0' numpy in a Python &gt;=3.10 interpreter; set ALPHAGENOME_API_KEY (and optionally ALPHAGENOME_PYTHON to that interpreter, ALPHAGENOME_TIMEOUT_MS). |
| Required environment | `ALPHAGENOME_API_KEY` (secret) |
| Effects | Sends variant coordinates and alleles through a Python subprocess to the Google DeepMind AlphaGenome API (Atlas lookups or live inference) using the caller's API key; no local files are written by the tools. The upstream README states the API is free for non-commercial use only, so commercial use is restricted and requires accepting DeepMind's AlphaGenome Terms of Use and rate limits. |

| Tool | Category | Summary |
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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `variants`
- **Keywords:** gwas locus, fine-mapping, variant prioritization, alphagenome, indel, quantile

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variants` | array of object | yes | Variants to score (1-100) |
| `chromosome` | string | no | Optional: locus chromosome, for the label |
| `start` | number | no | Optional: locus start, for the label |
| `end` | number | no | Optional: locus end, for the label |

??? note "Upstream description"

    Rank the variants of a locus by predicted effect (largest absolute quantile across modalities), to prioritize candidates for follow-up. For single-nucleotide variants only, atlas_lookup_variants or atlas_scan_region is faster and adds the AVI score.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alt`
- **Keywords:** alphagenome, regulatory variant, chromatin accessibility, histone marks, tf binding, splicing, hg38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | yes | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | no | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |

??? note "Upstream description"

    The predicted effects of a variant across every regulatory modality at once: accessibility, histone marks, TF binding, CAGE, RNA-seq, splice sites, polyadenylation and contact maps. Shows where the predicted effect concentrates; it does not label the variant.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alt`
- **Keywords:** alphagenome, avi score, variant effect, atlas, snv, indel, hg38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | yes | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | no | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |
| `source` | string | no | Optional: where the answer comes from (default: auto). auto = the precomputed AlphaGenome Atlas for single-nucleotide substitutions, live inference for everything else (indels, multi-nucleotide variants); falls back to live only when the Atlas does not hold the variant. atlas = Atlas only, errors instead of falling back. live = always run the model. The result always states which source answered. One of `auto`, `atlas`, `live` |
| `scorers` | array of string | no | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |

??? note "Upstream description"

    Predicted effect size of a variant across modalities, for prioritization. The tool name is kept for compatibility: it does NOT classify a variant as pathogenic or benign, and `classification` is always null.

    Returns the strongest effect per scorer (score, calibrated quantile, where it was seen), the largest absolute quantile, and, for a single-nucleotide variant answered from the Atlas, the AlphaGenome Variant Impact (AVI) score. `avi_score` is null on the live path, because the AVI score is served by the Atlas only.

    Source: a single-nucleotide variant is answered from the precomputed AlphaGenome Atlas; an indel or multi-nucleotide variant runs live inference (score_variant). Both return the same scorers in the same shape. Chosen automatically, overridable with `source`, and always stated in the result.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "How large is the predicted effect of chr19:44908684 T&gt;C?"

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alt`
- **Keywords:** alphagenome, atlas, snv, avi score, precomputed scores, hg38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference base (one of A, C, G, T). Must match hg38 at this position. |
| `alt` | string | yes | Alternate base (one of A, C, G, T) |
| `scorers` | array of string | no | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |
| `top_n` | number | no | Rows to return (default: 25, max: 100) |

??? note "Upstream description"

    Look up the precomputed AlphaGenome scores of one single-nucleotide variant. No model call, so it answers in seconds.

    Returns, per scorer, the strongest tracks for the variant ranked by absolute score, each with its calibrated quantile, gene, tissue or cell type, and assay. Default scorers: AVI_SCORE plus one per modality (RNA_SEQ, CAGE, DNASE, CHIP_HISTONE, CHIP_TF, SPLICE_SITES).

    If the reference base does not match hg38, the Atlas says which base it expected and that message is returned as a validation error.

    The Atlas holds precomputed AlphaGenome scores for single-nucleotide substitutions on the human reference genome (hg38, chr1-22, chrX, chrY). Indels and multi-nucleotide variants are not in it; use predict_variant_effect for those.

    The response is a summary, never a full score matrix: ranked rows only, capped at top_n (default 25, max 100) and at 40,000 characters.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Look up chr19:44908684 T&gt;C in the AlphaGenome Atlas"

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `variants`
- **Keywords:** alphagenome, atlas, batch snv, avi score, variant ranking, hg38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variants` | array of object | yes | Single-nucleotide variants to look up (1-500) |
| `scorers` | array of string | no | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |
| `top_n` | number | no | Rows to return (default: 25, max: 100) |

??? note "Upstream description"

    Look up precomputed AlphaGenome scores for up to 500 single-nucleotide variants in one call and rank them.

    Returns one row per variant (its strongest score and where it was seen), ranked. Default scorer: AVI_SCORE (AlphaGenome Variant Impact), one number per variant. With several scorers the ranking uses the largest absolute quantile. Variants the Atlas does not hold and variants it rejects (for example a reference base that does not match hg38) are listed separately with the reason; they do not fail the call.

    The Atlas holds precomputed AlphaGenome scores for single-nucleotide substitutions on the human reference genome (hg38, chr1-22, chrX, chrY). Indels and multi-nucleotide variants are not in it; use predict_variant_effect for those.

    The response is a summary, never a full score matrix: ranked rows only, capped at top_n (default 25, max 100) and at 40,000 characters.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Rank these 200 GWAS SNPs by their Atlas scores"

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `start`, `end`
- **Keywords:** alphagenome, saturation mutagenesis, region scan, atlas, avi score, hg38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `start` | number | yes | Start position (1-based, hg38) |
| `end` | number | yes | End position (greater than start; at most start + 10,000, or start + 50,000 with allow_large_region) |
| `allow_large_region` | boolean | no | Set to true to scan more than 10,000 bp (up to 50,000 bp). Slower, and the result may be incomplete (default: false) |
| `scorers` | array of string | no | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |
| `top_n` | number | no | Rows to return (default: 25, max: 100) |

??? note "Upstream description"

    Scan a genomic region in the AlphaGenome Atlas: every possible single-nucleotide substitution in the interval, ranked. Answers "which positions in this region matter most?" without running the model.

    Region width: at most 10,000 bp. Up to 50,000 bp only with allow_large_region=true. A scan is one API request per 32 bp under a requests-per-minute quota, so a large scan can take minutes. If the quota or the time limit stops a scan early, the partial result is returned, marked "Incomplete", with the range that was really scanned; the ranking then covers that range only.

    Default scorer: AVI_SCORE (about 7 seconds for 2,000 bp). Multi-track scorers are much slower (2,000 bp: DNASE 17 s, CHIP_TF 86 s). Scorers with one row per gene or junction (RNA_SEQ, SPLICE_JUNCTIONS, ...) cannot be used for a scan: scan with AVI_SCORE, then use atlas_lookup_variant on the top variants.

    The Atlas holds precomputed AlphaGenome scores for single-nucleotide substitutions on the human reference genome (hg38, chr1-22, chrX, chrY). Indels and multi-nucleotide variants are not in it; use predict_variant_effect for those.

    The response is a summary, never a full score matrix: ranked rows only, capped at top_n (default 25, max 100) and at 40,000 characters.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Scan chr17:49209289-49211289 and show the 10 substitutions with the largest predicted effect"

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `variants`, `modality`
- **Keywords:** alphagenome, variant ranking, expression, splicing, tf binding, chromatin

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variants` | array of object | yes | Variants to score (1-100) |
| `modality` | string | yes | One of `expression`, `splicing`, `tf_binding`, `chromatin` |

??? note "Upstream description"

    Rank several variants by predicted effect within one modality: expression (RNA_SEQ, CAGE), splicing (SPLICE_SITES, SPLICE_SITE_USAGE), tf_binding (CHIP_TF) or chromatin (DNASE, ATAC).

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `variants`
- **Keywords:** alphagenome, variant filter, quantile threshold, avi score, prioritization

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variants` | array of object | yes | Variants to score (1-100) |
| `threshold` | number | no | Smallest absolute quantile to keep, between 0 and 1 (default: 0.99) |
| `source` | string | no | Optional: where the answer comes from (default: auto). auto = the precomputed AlphaGenome Atlas for single-nucleotide substitutions, live inference for everything else (indels, multi-nucleotide variants); falls back to live only when the Atlas does not hold the variant. atlas = Atlas only, errors instead of falling back. live = always run the model. The result always states which source answered. One of `auto`, `atlas`, `live` |
| `scorers` | array of string | no | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |

??? note "Upstream description"

    Keep the variants whose predicted effect reaches a threshold, ranked. The tool name is kept for compatibility: it filters on predicted effect size, NOT on pathogenicity, and classifies nothing.

    `threshold` is the smallest absolute calibrated quantile (0 to 1) a variant must reach to be kept (default 0.99): the AVI score's quantile for variants answered from the Atlas, the largest quantile across modalities for live inference. Variants are routed per variant and reported per source; groups from different sources are not comparable.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Which of these variants have a predicted effect above the 99.9th percentile?"

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `variants`, `scoring_metric`
- **Keywords:** alphagenome, batch scoring, variant ranking, atlas, rna-seq, splicing

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variants` | array of object | yes | Variants to score (1-100) |
| `scoring_metric` | string | yes | What to rank by One of `rna_seq`, `splice`, `regulatory_impact`, `combined` |
| `top_n` | number | no | Variants to return per group (default: 10, max: 100) |
| `source` | string | no | Optional: where the answer comes from (default: auto). auto = the precomputed AlphaGenome Atlas for single-nucleotide substitutions, live inference for everything else (indels, multi-nucleotide variants); falls back to live only when the Atlas does not hold the variant. atlas = Atlas only, errors instead of falling back. live = always run the model. The result always states which source answered. One of `auto`, `atlas`, `live` |
| `scorers` | array of string | no | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |

??? note "Upstream description"

    Score up to 100 variants and rank them by predicted effect.

    Each variant is routed on its own: single-nucleotide variants to the Atlas, the rest to live inference. A mixed batch comes back as two separately ranked groups, because the Atlas group is ranked by the AVI score and live inference has no AVI score; the two must not be compared. The result reports how many variants came from each source and how many fell back.

    Scoring metric (used when `scorers` is not given): rna_seq = RNA_SEQ, splice = SPLICE_SITES, regulatory_impact and combined = AVI_SCORE from the Atlas and every modality (ranked by the largest absolute quantile) from live inference.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Score these 50 variants and show me the top 10 by predicted effect"

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `variants`, `tissues`
- **Keywords:** alphagenome, tissue-specific, variant ranking, uberon, tissues

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variants` | array of object | yes | Variants to score (1-100) |
| `tissues` | array of string | yes | Tissue names or ontology CURIEs (1-10) |

??? note "Upstream description"

    Rank several variants by predicted effect within each of several tissues: one ranking per tissue.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alts`
- **Keywords:** alphagenome, alternate alleles, allele comparison, hg38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference allele (A, C, G, T; more than one base for an indel) |
| `alts` | array of string | yes | Alternate alleles (1-20) |

??? note "Upstream description"

    Rank the alternate alleles of one position by predicted effect.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `protective_variant`, `risk_variant`
- **Keywords:** alphagenome, protective allele, risk allele, variant comparison

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `protective_variant` | object | yes |  |
| `risk_variant` | object | yes |  |

??? note "Upstream description"

    Two variants side by side, labelled as the caller names them. The tool compares predicted effect sizes per modality; it does not judge which allele is protective or a risk.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `variant1`, `variant2`
- **Keywords:** alphagenome, variant comparison, rsid, modality, quantile

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant1` | object | yes |  |
| `variant2` | object | yes |  |

??? note "Upstream description"

    Two variants side by side: the strongest predicted effect of each modality for both, and which of the two has the larger absolute quantile per scorer. A comparison of predicted effect sizes, not of severity.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Compare APOE rs429358 and rs7412"

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `variants`
- **Keywords:** alphagenome, gene symbol, variant ranking, splicing, expression

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variants` | array of object | yes | Variants to score (1-100) |
| `gene_name` | string | no | Optional: gene symbol (e.g., APOE) |

??? note "Upstream description"

    Rank several variants by their predicted effect on one gene. With gene_name, the gene-level scorers (RNA_SEQ, SPLICE_SITES) are restricted to that gene.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alt`
- **Keywords:** alphagenome, variant interpretation, avi score, natural language summary, hg38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | yes | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | no | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |
| `source` | string | no | Optional: where the answer comes from (default: auto). auto = the precomputed AlphaGenome Atlas for single-nucleotide substitutions, live inference for everything else (indels, multi-nucleotide variants); falls back to live only when the Atlas does not hold the variant. atlas = Atlas only, errors instead of falling back. live = always run the model. The result always states which source answered. One of `auto`, `atlas`, `live` |
| `scorers` | array of string | no | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |

??? note "Upstream description"

    Plain sentences that restate a variant's predicted effects: the AVI score and its largest contributions (from the Atlas), then the strongest effect of each modality ordered by absolute quantile, with the direction for signed scorers. Descriptive only: the sentences restate returned numbers and make no statement about pathogenicity.

    Source: a single-nucleotide variant is answered from the precomputed AlphaGenome Atlas; an indel or multi-nucleotide variant runs live inference (score_variant). Both return the same scorers in the same shape. Chosen automatically, overridable with `source`, and always stated in the result.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Explain the predicted effect of chr17:49210289 C&gt;T in plain language"

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alt`
- **Keywords:** alphagenome, variant report, avi score, feature importance, hg38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | yes | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | no | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |
| `source` | string | no | Optional: where the answer comes from (default: auto). auto = the precomputed AlphaGenome Atlas for single-nucleotide substitutions, live inference for everything else (indels, multi-nucleotide variants); falls back to live only when the Atlas does not hold the variant. atlas = Atlas only, errors instead of falling back. live = always run the model. The result always states which source answered. One of `auto`, `atlas`, `live` |
| `scorers` | array of string | no | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |

??? note "Upstream description"

    A fuller report of one variant's predicted molecular effects: more rows per scorer than predict_variant_effect and, from the Atlas, the AVI score with its feature attributions (AVI_SCORE_FEATURE_IMPORTANCE). It is a research summary, not a clinical report: it contains no pathogenicity classification and no recommendation.

    Source: a single-nucleotide variant is answered from the precomputed AlphaGenome Atlas; an indel or multi-nucleotide variant runs live inference (score_variant). Both return the same scorers in the same shape. Chosen automatically, overridable with `source`, and always stated in the result.

    The response is a summary, never a full score matrix: ranked rows only, capped at top_n (default 25, max 100) and at 40,000 characters.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Generate a report for chr19:44908684 T&gt;C"

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alt`
- **Keywords:** alphagenome, allele-specific expression, log fold change, rna-seq, eqtl-like

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | yes | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | no | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |

??? note "Upstream description"

    Predicted expression with the alternate allele against the reference allele: the RNA_SEQ scorer is that log fold change per gene and tissue, and RNA_SEQ_ACTIVE gives the expression level alongside it.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alt`
- **Keywords:** alphagenome, atac-seq, dnase-seq, chromatin accessibility, regulatory variant

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | yes | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | no | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |

??? note "Upstream description"

    Predicted chromatin accessibility effects of a variant (ATAC-seq and DNase-seq).

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alt`
- **Keywords:** alphagenome, gene expression, rna-seq, cage, log fold change

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | yes | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | no | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |

??? note "Upstream description"

    Predicted gene expression effects of a variant: RNA-seq (log fold change per gene and tissue) and CAGE.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alt`
- **Keywords:** alphagenome, splicing, splice junction, splice site usage, sqtl

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | yes | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | no | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |

??? note "Upstream description"

    Predicted splicing effects of a variant: splice sites, splice site usage and splice junctions, with the gene and junction of each.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alt`
- **Keywords:** alphagenome, tf chip-seq, transcription factor, binding site, regulatory variant

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | yes | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | no | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |

??? note "Upstream description"

    Predicted transcription factor binding effects of a variant (TF ChIP-seq), with the factor and cell type of each.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alt`
- **Keywords:** alphagenome, tissue-specific, brain, liver, heart, variant effect

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | yes | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissues` | array of string | no | Tissue names (brain, neuron, blood, liver, heart, lung, kidney) or ontology CURIEs. Default: brain, liver, heart |

??? note "Upstream description"

    The strongest predicted effect of a variant in each of several tissues: the same scores, filtered to the tracks of one tissue at a time.

    Runs live inference (score_variant with the SDK's recommended variant scorers); works for single-nucleotide variants, indels and multi-nucleotide variants. The result states `source: live`.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Compare the predicted effect of chr19:44908684 T&gt;C in brain, liver and heart"

Call it through DNAKit once the server is installed:

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

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `ref`, `alt`
- **Keywords:** alphagenome, variant effect prediction, avi score, regulatory, splicing, hg38

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome (chr1-chr22, chrX, chrY) |
| `position` | number | yes | Genomic position (1-based, hg38) |
| `ref` | string | yes | Reference allele (A, C, G, T; more than one base for an indel) |
| `alt` | string | yes | Alternate allele (A, C, G, T; more than one base for an indel) |
| `tissue_type` | string | no | Optional: keep only the tracks of one tissue or cell type. A name (brain, neuron, blood, liver, heart, lung, kidney) or an ontology CURIE (e.g., UBERON:0000955, CL:0000540). Default: all tissues. |
| `output_types` | array of string | no | Optional: modalities to report (default: one scorer per modality) |
| `source` | string | no | Optional: where the answer comes from (default: auto). auto = the precomputed AlphaGenome Atlas for single-nucleotide substitutions, live inference for everything else (indels, multi-nucleotide variants); falls back to live only when the Atlas does not hold the variant. atlas = Atlas only, errors instead of falling back. live = always run the model. The result always states which source answered. One of `auto`, `atlas`, `live` |
| `scorers` | array of string | no | Optional: scorer names to use instead of the defaults. Names come from atlas_list_scorers and are the same for both sources, except the AVI scorers, which the Atlas alone serves. |

??? note "Upstream description"

    Predicted regulatory effect of a genetic variant, per modality: the strongest tracks of each scorer (expression, transcription start, chromatin accessibility, histone marks, transcription factor binding, splicing), each with its score, calibrated quantile, gene, tissue or cell type, and assay. From the Atlas the AlphaGenome Variant Impact (AVI) score is included.

    Source: a single-nucleotide variant is answered from the precomputed AlphaGenome Atlas; an indel or multi-nucleotide variant runs live inference (score_variant). Both return the same scorers in the same shape. Chosen automatically, overridable with `source`, and always stated in the result.

    The response is a summary, never a full score matrix: ranked rows only, capped at top_n (default 25, max 100) and at 40,000 characters.

    Results are AlphaGenome model predictions for research prioritization, not clinical classifications: scores and calibrated quantiles are reported as returned, and no pathogenic/benign call is made.

    Example: "Analyze chr19:44908684 T&gt;C with AlphaGenome"

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'alphagenome-jolab',
    'predict_variant_effect',
    {'chromosome': '...', 'position': 0.0, 'ref': '...', 'alt': '...'},
)
```
