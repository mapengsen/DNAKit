# ToolUniverse: prediction_model (`tooluniverse`)

DNAKit's MCP tool index keeps **845** DNA/gene tools of this server (2927 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

This page lists the 17 tools in category `prediction_model`; other categories: [`sequence`](tooluniverse-sequence.md) (44), [`genome_annotation`](tooluniverse-genome_annotation.md) (253), [`variant`](tooluniverse-variant.md) (155), [`association`](tooluniverse-association.md) (138), [`expression`](tooluniverse-expression.md) (56), [`regulation_epigenome`](tooluniverse-regulation_epigenome.md) (54), [`sequencing_data`](tooluniverse-sequencing_data.md) (91), [`molecular_tools`](tooluniverse-molecular_tools.md) (23), [`breeding_population`](tooluniverse-breeding_population.md) (14).

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
| [`AlphaGenome_atlas_lookup_variant`](#AlphaGenome_atlas_lookup_variant) | `prediction_model` | Looks up precomputed AlphaGenome Atlas effect scores for a single human SNV. |
| [`AlphaGenome_atlas_scan_interval`](#AlphaGenome_atlas_scan_interval) | `prediction_model` | Returns precomputed AlphaGenome Atlas scores for all possible SNVs in a region up to 10 kb. |
| [`AlphaGenome_predict_interval`](#AlphaGenome_predict_interval) | `prediction_model` | Predicts multimodal genomic tracks (expression, accessibility, etc.) for a genomic interval with AlphaGenome. |
| [`AlphaGenome_predict_sequence`](#AlphaGenome_predict_sequence) | `prediction_model` | Predicts multimodal genomic tracks for a raw DNA sequence with AlphaGenome. |
| [`AlphaGenome_predict_variant`](#AlphaGenome_predict_variant) | `prediction_model` | Predicts full multimodal tracks for reference and alternate alleles of a variant with AlphaGenome. |
| [`AlphaGenome_score_interval`](#AlphaGenome_score_interval) | `prediction_model` | Scores baseline per-gene, per-track activity for a genomic interval with AlphaGenome. |
| [`AlphaGenome_score_ism_variants`](#AlphaGenome_score_ism_variants) | `prediction_model` | Scores all single-nucleotide substitutions in a window up to 500 bp with AlphaGenome and ranks them. |
| [`AlphaGenome_score_variant`](#AlphaGenome_score_variant) | `prediction_model` | Scores a variant's ref-vs-alt effect per track with AlphaGenome and returns top effects. |
| [`ESM2_score_missense_variant`](#ESM2_score_missense_variant) | `prediction_model` | Scores a missense variant with the ESM-2 masked-marginal log-likelihood ratio. |
| [`ESM_explain_variant_mechanism`](#ESM_explain_variant_mechanism) | `prediction_model` | Scores a missense variant with ESMC-6B SAE feature disruption and summarizes affected feature categories. |
| [`ESM_score_variant_sae_batch`](#ESM_score_variant_sae_batch) | `prediction_model` | Scores many missense variants of one protein with ESMC-6B SAE disruption. |
| [`ESM_score_variant_sae_disruption`](#ESM_score_variant_sae_disruption) | `prediction_model` | Scores a missense variant by ESMC-6B SAE feature activation changes around the site. |
| [`Evo2_score_variant`](#Evo2_score_variant) | `prediction_model` | Scores a DNA variant zero-shot with Evo 2 as the alt-minus-ref log-likelihood. |
| [`NvidiaNIM_evo2`](#NvidiaNIM_evo2) | `prediction_model` | Generates DNA sequence continuations with the Evo 2 genome language model (40B or 7B) via NVIDIA NIM. |
| [`SpliceAI_get_max_delta`](#SpliceAI_get_max_delta) | `prediction_model` | Returns the maximum SpliceAI delta score across splice effects for a variant with a benign-to-high pathogenicity interpretation. |
| [`SpliceAI_predict_pangolin`](#SpliceAI_predict_pangolin) | `prediction_model` | Predicts splice-altering effects of a variant per transcript with the Pangolin deep learning model. |
| [`SpliceAI_predict_splice`](#SpliceAI_predict_splice) | `prediction_model` | Predicts splice-altering effects of a variant with SpliceAI, returning acceptor/donor gain and loss delta scores. |

## `AlphaGenome_atlas_lookup_variant` {#AlphaGenome_atlas_lookup_variant}

Looks up precomputed AlphaGenome Atlas effect scores for a single human SNV.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `reference_bases`, `alternate_bases`
- **Keywords:** snv, variant effect, avi score, precomputed, chromosome position
- **Needs API keys:** `ALPHA_GENOME_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome, e.g. 'chr22'. |
| `position` | integer | yes | 1-based variant position. |
| `reference_bases` | string | yes | Reference allele, e.g. 'A'. |
| `alternate_bases` | string | yes | Alternate allele, e.g. 'C'. |
| `scorers` | array \| null | no | Which precomputed Atlas scorers to return, e.g. \['AVI_SCORE', 'RNA_SEQ', 'ATAC', 'CHIP_TF'\]. Defaults to \['AVI_SCORE'\], the unified variant impact score. |
| `top_n` | integer \| null | no | Number of top \|effect\| tracks to return per scorer (default 20). |

??? note "Upstream description"

    Look up precomputed variant effect scores for a single-nucleotide substitution from AlphaGenome Atlas (Google DeepMind, released Sep 2026): a precomputed catalogue covering all ~9 billion possible human SNVs, including the unified AlphaGenome Variant Impact (AVI) score that combines AlphaGenome and AlphaMissense predictions into one coding+non-coding pathogenicity-style number. This is a database lookup (fast, high query rate), not a live model run -- use AlphaGenome_score_variant instead for indels or live per-track scoring. Hosted API via the alphagenome SDK; free for non-commercial use (set ALPHA_GENOME_API_KEY).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AlphaGenome_atlas_lookup_variant',
     'arguments': {'chromosome': '...',
                   'position': 0,
                   'reference_bases': '...',
                   'alternate_bases': '...'}},
)
```

## `AlphaGenome_atlas_scan_interval` {#AlphaGenome_atlas_scan_interval}

Returns precomputed AlphaGenome Atlas scores for all possible SNVs in a region up to 10 kb.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `start`, `end`
- **Keywords:** saturation scan, snv, promoter, variant ranking, genomic region
- **Needs API keys:** `ALPHA_GENOME_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome, e.g. 'chr22'. |
| `start` | integer | yes | Interval start (0-based). Window (end - start) must be &lt;= 10000 bp. |
| `end` | integer | yes | Interval end. Window (end - start) must be &lt;= 10000 bp. |
| `scorers` | array \| null | no | Which precomputed Atlas scorers to return (see AlphaGenome_atlas_list_scorers), e.g. \['AVI_SCORE'\]. Defaults to \['AVI_SCORE'\]. |
| `top_n` | integer \| null | no | Number of top \|effect\| entries to return (default 20). |

??? note "Upstream description"

    Look up precomputed AlphaGenome Atlas scores for EVERY possible single-nucleotide substitution in a genomic region (&lt;=10,000 bp), e.g. to rank all variants in a promoter by AVI score at once. This is a precomputed database read (fast, high query rate), not a live model run -- for a single known variant use AlphaGenome_atlas_lookup_variant instead. Hosted API via the alphagenome SDK; free for non-commercial use (set ALPHA_GENOME_API_KEY).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AlphaGenome_atlas_scan_interval',
     'arguments': {'chromosome': '...', 'start': 0, 'end': 0}},
)
```

## `AlphaGenome_predict_interval` {#AlphaGenome_predict_interval}

Predicts multimodal genomic tracks (expression, accessibility, etc.) for a genomic interval with AlphaGenome.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `start`, `end`
- **Keywords:** genomic tracks, rna-seq, atac, chip, interval prediction
- **Needs API keys:** `ALPHA_GENOME_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome, e.g. 'chr19'. |
| `start` | integer | yes | Interval start (0-based). |
| `end` | integer | yes | Interval end. |
| `output_types` | array \| null | no | Modalities to predict (default \['RNA_SEQ'\]): RNA_SEQ, ATAC, DNASE, CAGE, CHIP_HISTONE, CHIP_TF, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, CONTACT_MAPS, PROCAP. |
| `ontology_terms` | array \| null | no | Optional tissue/cell ontology terms (e.g. \['UBERON:0001114'\] = liver) to restrict tracks. |
| `organism` | string \| null | no | 'human' (default) or 'mouse'. |
| `sequence_length` | string \| null | no | Context window: 16KB, 100KB, 500KB, or 1MB (default). |

??? note "Upstream description"

    Predict multimodal genomic tracks for a genomic interval with DeepMind AlphaGenome (single DNA-sequence model; up to 1 Mb at single-base resolution). Returns a compact per-modality summary (track counts and shapes) for the requested outputs. Hosted API via the alphagenome SDK; free for non-commercial use (set ALPHA_GENOME_API_KEY).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AlphaGenome_predict_interval',
     'arguments': {'chromosome': '...', 'start': 0, 'end': 0}},
)
```

## `AlphaGenome_predict_sequence` {#AlphaGenome_predict_sequence}

Predicts multimodal genomic tracks for a raw DNA sequence with AlphaGenome.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequence`
- **Keywords:** dna sequence, synthetic construct, track prediction
- **Needs API keys:** `ALPHA_GENOME_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Raw DNA sequence (A/C/G/T/N). Length must equal a supported context window (16384/131072/524288/1048576 bp). |
| `output_types` | array \| null | no | Modalities to predict (default \['RNA_SEQ'\]): RNA_SEQ, ATAC, DNASE, CAGE, CHIP_HISTONE, CHIP_TF, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, CONTACT_MAPS, PROCAP. |
| `ontology_terms` | array \| null | no | Optional tissue/cell ontology terms (e.g. \['UBERON:0001114'\] = liver) to restrict tracks. Use AlphaGenome_output_metadata to discover valid values. |
| `organism` | string \| null | no | 'human' (default) or 'mouse'. |

??? note "Upstream description"

    Predict multimodal genomic tracks for a raw DNA sequence (no genome coordinates needed) with DeepMind AlphaGenome. Useful for synthetic constructs or sequences not anchored to a reference genome. The sequence length must match one of the model's supported context windows (16384, 131072, 524288, or 1048576 bp; pad with 'N' as needed). Hosted API via the alphagenome SDK; free for non-commercial use (set ALPHA_GENOME_API_KEY).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AlphaGenome_predict_sequence',
     'arguments': {'sequence': '...'}},
)
```

## `AlphaGenome_predict_variant` {#AlphaGenome_predict_variant}

Predicts full multimodal tracks for reference and alternate alleles of a variant with AlphaGenome.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `reference_bases`, `alternate_bases`
- **Keywords:** variant effect, ref alt, genomic tracks, regulatory variant
- **Needs API keys:** `ALPHA_GENOME_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome, e.g. 'chr22'. |
| `position` | integer | yes | 1-based variant position. |
| `reference_bases` | string | yes | Reference allele, e.g. 'A'. |
| `alternate_bases` | string | yes | Alternate allele, e.g. 'C'. |
| `output_types` | array \| null | no | Modalities to predict (default \['RNA_SEQ'\]): RNA_SEQ, ATAC, DNASE, CAGE, CHIP_HISTONE, CHIP_TF, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, CONTACT_MAPS, PROCAP. |
| `ontology_terms` | array \| null | no | Optional tissue/cell ontology terms (e.g. \['UBERON:0001114'\] = liver) to restrict tracks. Use AlphaGenome_output_metadata to discover valid values. |
| `organism` | string \| null | no | 'human' (default) or 'mouse'. |
| `sequence_length` | string \| null | no | Context window: 16KB, 100KB, 500KB, or 1MB (default). |

??? note "Upstream description"

    Predict full multimodal tracks for BOTH the reference and alternate allele of a variant with DeepMind AlphaGenome (distinct from AlphaGenome_score_variant, which reduces this to a single recommended per-gene effect score). Returns a compact per-modality summary (shape/track counts) for each allele -- use this when you need the raw predicted tracks for both alleles rather than a scored comparison. Hosted API via the alphagenome SDK; free for non-commercial use (set ALPHA_GENOME_API_KEY).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AlphaGenome_predict_variant',
     'arguments': {'chromosome': '...',
                   'position': 0,
                   'reference_bases': '...',
                   'alternate_bases': '...'}},
)
```

## `AlphaGenome_score_interval` {#AlphaGenome_score_interval}

Scores baseline per-gene, per-track activity for a genomic interval with AlphaGenome.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `start`, `end`
- **Keywords:** gene expression prediction, interval score, gene mask
- **Needs API keys:** `ALPHA_GENOME_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome, e.g. 'chr22'. |
| `start` | integer | yes | Interval start (0-based). |
| `end` | integer | yes | Interval end. |
| `organism` | string \| null | no | 'human' (default) or 'mouse'. |
| `sequence_length` | string \| null | no | Context window: 16KB, 100KB, 500KB, or 1MB (default). |
| `top_n` | integer \| null | no | Number of top \|effect\| entries to return (default 20). |

??? note "Upstream description"

    Score the default gene-mask activity for a genomic interval with DeepMind AlphaGenome, without reference to a specific variant (e.g. baseline predicted expression level per nearby gene). Returns the top per-gene, per-track scores sorted by \|effect\|. For variant-effect scoring (ref vs alt), use AlphaGenome_score_variant instead. Hosted API via the alphagenome SDK; free for non-commercial use (set ALPHA_GENOME_API_KEY).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AlphaGenome_score_interval',
     'arguments': {'chromosome': '...', 'start': 0, 'end': 0}},
)
```

## `AlphaGenome_score_ism_variants` {#AlphaGenome_score_ism_variants}

Scores all single-nucleotide substitutions in a window up to 500 bp with AlphaGenome and ranks them.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `start`, `end`
- **Keywords:** in silico mutagenesis, ism, saturation mutagenesis, regulatory bases
- **Needs API keys:** `ALPHA_GENOME_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome, e.g. 'chr22'. |
| `start` | integer | yes | ISM window start (0-based). Window (end - start) must be &lt;= 500 bp. |
| `end` | integer | yes | ISM window end. Window (end - start) must be &lt;= 500 bp. |
| `output_type` | string \| null | no | Modality/scorer to rank by, i.e. a key of alphagenome.models.variant_scorers.RECOMMENDED_VARIANT_SCORERS. RNA_SEQ (default), ATAC, DNASE, CAGE, CHIP_HISTONE, CHIP_TF, PROCAP, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, POLYADENYLATION, CONTACT_MAPS, or an _ACTIVE variant of ATAC/DNASE/CAGE/CHIP_TF/CHIP_HISTONE/PROCAP/RNA_SEQ. |
| `organism` | string \| null | no | 'human' (default) or 'mouse'. |
| `sequence_length` | string \| null | no | Context window: 16KB, 100KB, 500KB, or 1MB (default). |
| `top_n` | integer \| null | no | Number of top \|effect\| entries to return (default 20). |

??? note "Upstream description"

    In-silico saturation mutagenesis with DeepMind AlphaGenome: scores every possible single-nucleotide substitution across a short window (&lt;=500 bp) and ranks candidates by their single most extreme per-gene, per-track effect -- a 'which bases matter most' scan for finding causal regulatory positions (e.g. within a promoter or enhancer). Each live call costs ~3x the window width in underlying model evaluations, so the window is capped at 500 bp; for a single known variant use AlphaGenome_score_variant instead. Hosted API via the alphagenome SDK; free for non-commercial use (set ALPHA_GENOME_API_KEY).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AlphaGenome_score_ism_variants',
     'arguments': {'chromosome': '...', 'start': 0, 'end': 0}},
)
```

## `AlphaGenome_score_variant` {#AlphaGenome_score_variant}

Scores a variant's ref-vs-alt effect per track with AlphaGenome and returns top effects.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `chromosome`, `position`, `reference_bases`, `alternate_bases`
- **Keywords:** variant effect, regulatory variant, eqtl prediction, sequence model
- **Needs API keys:** `ALPHA_GENOME_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `chromosome` | string | yes | Chromosome, e.g. 'chr22'. |
| `position` | integer | yes | 1-based variant position. |
| `reference_bases` | string | yes | Reference allele, e.g. 'A'. |
| `alternate_bases` | string | yes | Alternate allele, e.g. 'C'. |
| `output_type` | string \| null | no | Modality to score, i.e. a key of alphagenome.models.variant_scorers.RECOMMENDED_VARIANT_SCORERS. RNA_SEQ (default), ATAC, DNASE, CAGE, CHIP_HISTONE, CHIP_TF, PROCAP, SPLICE_SITES, SPLICE_SITE_USAGE, SPLICE_JUNCTIONS, POLYADENYLATION, CONTACT_MAPS, or an _ACTIVE variant of ATAC/DNASE/CAGE/CHIP_TF/CHIP_HISTONE/PROCAP/RNA_SEQ. |
| `organism` | string \| null | no | 'human' (default) or 'mouse'. |
| `sequence_length` | string \| null | no | Context window: 16KB, 100KB, 500KB, or 1MB (default). |
| `top_n` | integer \| null | no | Number of top \|effect\| tracks to return (default 20). |

??? note "Upstream description"

    Score a regulatory variant's effect with DeepMind AlphaGenome (Avsec, Nature 2026), the hosted successor to Enformer/Borzoi. Predicts ref vs alt over up to 1 Mb at single-base resolution and returns the recommended per-track effect scores (sorted by \|effect\|). Hosted API via the alphagenome SDK; free for non-commercial use (set ALPHA_GENOME_API_KEY).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'AlphaGenome_score_variant',
     'arguments': {'chromosome': '...',
                   'position': 0,
                   'reference_bases': '...',
                   'alternate_bases': '...'}},
)
```

## `ESM2_score_missense_variant` {#ESM2_score_missense_variant}

Scores a missense variant with the ESM-2 masked-marginal log-likelihood ratio.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequence`, `position`, `mutant`
- **Keywords:** esm-2, missense, llr, variant effect
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Wild-type protein sequence in 1-letter amino-acid code (whitespace is ignored). |
| `position` | integer | yes | 1-based residue position of the variant within the sequence. |
| `mutant` | string | yes | Mutant amino acid (single 1-letter code, one of the 20 standard residues). |
| `wild_type` | string \| null | no | Optional wild-type residue (1-letter). If given, it is validated against the residue at `position` to catch coordinate/sequence mismatches. |
| `model_id` | string \| null | no | ESM-2 model to use. Default 'facebook/esm2_t33_650M_UR50D'. Other sizes: esm2_t30_150M_UR50D, esm2_t12_35M_UR50D, esm2_t6_8M_UR50D. |
| `wait_for_model` | boolean \| null | no | If true, block until the model finishes loading on the HF servers instead of returning a retryable 'loading' status. |

??? note "Upstream description"

    Score a missense protein variant with ESM-2 masked-marginal log-likelihood ratio (Meier 2021) over HuggingFace's keyless hf-inference. Input wild-type 1-letter sequence + 1-based position + mutant residue; returns LLR = logP(mut) - logP(wt). Negative =&gt; mutant disfavored (candidate deleterious). No API key. For ESMC/SAE scoring with a key use ESM_score_sequence.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ESM2_score_missense_variant',
     'arguments': {'sequence': '...', 'position': 0, 'mutant': '...'}},
)
```

## `ESM_explain_variant_mechanism` {#ESM_explain_variant_mechanism}

Scores a missense variant with ESMC-6B SAE feature disruption and summarizes affected feature categories.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequence`, `position`, `ref_aa`, `alt_aa`
- **Keywords:** missense, variant mechanism, esmc, sae
- **Needs API keys:** `ESM_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Reference (wild-type) protein amino acid sequence in single-letter code. Up to ~2700 AA. |
| `position` | integer | yes | 1-indexed residue position of the variant |
| `ref_aa` | string | yes | Single-letter wild-type amino acid (must match sequence\[position-1\]) |
| `alt_aa` | string | yes | Single-letter substituted amino acid |
| `window` | integer | no | Residue window centered on the mutation for activation summation. Default `8` |
| `top_k_features` | integer | no | Number of top lost / top gained features to describe and include in the category summary. Default `5` |
| `include_descriptions` | boolean | no | If true (default), call ESM_describe_sae_feature on each top feature to get category labels. If false, return raw feature IDs without categories (cheaper). Default `True` |
| `model` | string | no | ESMC base model One of `esmc-6b-2024-12` Default `esmc-6b-2024-12` |
| `sae_model` | string | no | SAE codebook identifier Default `esmc-6b-2024-12_k64_codebook16384_layer60` |

??? note "Upstream description"

    One-call composite for variant mechanism: runs ESMC-6B SAE variant disruption + describe_sae_feature on each top affected feature + composes a category-level mechanism summary (e.g. 'Disrupted feature categories (lost): catalytic=2, ligand-binding=1'). Use when a variant-interpretation skill wants the mechanism in a single tool call rather than orchestrating disruption + N description calls. Set include_descriptions=false to skip the category labeling step (2 Forge calls only). Otherwise: 2 disruption calls + up to (2 * top_k_features) describe calls (most cached after first use). Prerequisites: pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52'; set ESM_API_KEY. Outputs governed by EvolutionaryScale Cambrian Inference License — non-commercial use only unless covered by a separate commercial agreement.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ESM_explain_variant_mechanism',
     'arguments': {'sequence': '...',
                   'position': 0,
                   'ref_aa': '...',
                   'alt_aa': '...'}},
)
```

## `ESM_score_variant_sae_batch` {#ESM_score_variant_sae_batch}

Scores many missense variants of one protein with ESMC-6B SAE disruption.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequence`, `variants`
- **Keywords:** missense, variant scoring, batch, esmc
- **Needs API keys:** `ESM_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Reference (wild-type) protein amino acid sequence in single-letter code. Up to ~2700 AA. |
| `variants` | array of object | yes | List of variants to score. Each variant: \{position: 1-indexed int, ref_aa: single letter matching sequence\[position-1\], alt_aa: single letter substitution\}. |
| `window` | integer | no | Residue window centered on each variant for activation summation. Default 8 (i.e. positions \[pos-8, pos+8\]). Default `8` |
| `top_k_features` | integer | no | Number of top features to return per variant (top-K lost + top-K gained). Default `10` |
| `model` | string | no | ESMC base model (must be the 6B model — SAEs are trained against it) One of `esmc-6b-2024-12` Default `esmc-6b-2024-12` |
| `sae_model` | string | no | SAE codebook identifier Default `esmc-6b-2024-12_k64_codebook16384_layer60` |

??? note "Upstream description"

    Score many missense variants against one reference protein using ESMC-6B SAE. More Forge-efficient than calling ESM_score_variant_sae_disruption N times: runs the reference SAE once and reuses it across every variant, giving N+1 Forge calls instead of 2N. Use for saturation mutagenesis at a single position (all 19 alternative residues), DMS-style sweeps, clinical variant panels, or any workflow scoring multiple variants on the same protein. Max 100 variants per call. Prerequisites: pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52' (SAE support is on an unmerged feature branch — PyPI esm 3.2.x does NOT include SAEConfig). Set ESM_API_KEY environment variable. Outputs governed by EvolutionaryScale Cambrian Inference License — non-commercial use only unless covered by a separate commercial agreement.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ESM_score_variant_sae_batch',
     'arguments': {'sequence': '...', 'variants': []}},
)
```

## `ESM_score_variant_sae_disruption` {#ESM_score_variant_sae_disruption}

Scores a missense variant by ESMC-6B SAE feature activation changes around the site.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequence`, `position`, `ref_aa`, `alt_aa`
- **Keywords:** missense, variant effect, esmc, sae
- **Needs API keys:** `ESM_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Reference protein sequence (canonical isoform). Single-letter codes, no gaps. The mutant sequence is built internally by substituting alt_aa at position. |
| `position` | integer | yes | 1-indexed mutation position. The amino acid at sequence\[position-1\] must equal ref_aa; otherwise the tool returns an error (the wrong sequence was supplied or the variant notation references a different isoform). |
| `ref_aa` | string | yes | Reference amino acid at the mutation position, single-letter code (e.g. 'R'). Tool validates this matches the supplied sequence. |
| `alt_aa` | string | yes | Mutant amino acid, single-letter code (e.g. 'H' for R175H). |
| `window` | integer | no | Residue window radius around the mutation. Per-feature activations are summed across this window before computing the delta. Default 8 = +/- 8 residues = 17-residue window. Default `8` |
| `top_k_features` | integer | no | Number of top LOST and top GAINED features to return. Default 10. Default `10` |
| `model` | string | no | ESMC backbone (default esmc-6b-2024-12). One of `esmc-6b-2024-12` Default `esmc-6b-2024-12` |
| `sae_model` | string | no | SAE checkpoint (default layer-60 6B SAE). Default `esmc-6b-2024-12_k64_codebook16384_layer60` |

??? note "Upstream description"

    Composite SAE-based variant scoring. Given a protein sequence and a missense variant (position + ref_aa + alt_aa), runs ESMC-6B SAE on both reference and mutant sequences, computes per-feature activation deltas summed over a window centered on the mutation site, and returns ranked top-K features LOST and GAINED in the mutant. This is the convenience layer over ESM_get_sae_features — use this for one-shot variant interpretation; use ESM_get_sae_features directly only when you need raw per-residue features for non-variant analyses. Validates ref_aa matches the position in the supplied sequence (returns clear error if not, so you can detect transcript / isoform mismatches). Cost: 2 Forge API credits (1 ref + 1 mut). Latency: ~3-6s for typical-length human proteins. Prerequisites: same as ESM_get_sae_features (pip install 'esm @ git+https://github.com/evolutionaryscale/esm@ee891c52', ESM_API_KEY env var). License: SAE outputs governed by Cambrian Inference License — non-commercial / academic only.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'ESM_score_variant_sae_disruption',
     'arguments': {'sequence': '...',
                   'position': 0,
                   'ref_aa': '...',
                   'alt_aa': '...'}},
)
```

## `Evo2_score_variant` {#Evo2_score_variant}

Scores a DNA variant zero-shot with Evo 2 as the alt-minus-ref log-likelihood.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** none
- **Keywords:** evo 2, genome language model, zero-shot, variant effect
- **Needs API keys:** `NVIDIA_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `ref_sequence` | string \| null | no | Reference DNA window (A/C/G/T/N). Use with alt_sequence (same length). |
| `alt_sequence` | string \| null | no | Alternate DNA window, same length/centering as ref_sequence. |
| `sequence` | string \| null | no | Reference DNA window for point-substitution mode (use with position + alternate). |
| `position` | integer \| null | no | 1-based position of the substituted base within `sequence`. |
| `reference` | string \| null | no | Optional reference base (single letter) at `position`, validated against `sequence`. |
| `alternate` | string \| null | no | Alternate base (single letter) substituted at `position`. |
| `model` | string | no | Hosted Evo 2 model size: evo2-40b (default, most accurate) or evo2-7b (smaller, faster). One of `evo2-40b`, `evo2-7b` Default `evo2-40b` |

??? note "Upstream description"

    Zero-shot variant-effect scoring with NVIDIA-hosted Evo 2 (Arc Institute genome foundation model). Runs the model's forward pass on the reference and alternate DNA windows and returns the delta log-likelihood (alt - ref). Negative = variant disfavored by the genome model (candidate deleterious). Two hosted sizes are selectable via `model`: evo2-40b (default, most accurate) and evo2-7b (smaller, faster). Hosted NIM API; requires NVIDIA_API_KEY. Complements NvidiaNIM_evo2 (which only generates sequences).

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'Evo2_score_variant', 'arguments': {}},
)
```

## `NvidiaNIM_evo2` {#NvidiaNIM_evo2}

Generates DNA sequence continuations with the Evo 2 genome language model (40B or 7B) via NVIDIA NIM.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequence`
- **Keywords:** evo 2, dna language model, sequence generation, logits
- **Needs API keys:** `NVIDIA_API_KEY`
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes | Input DNA sequence (A, C, T, G characters only) |
| `model` | string | no | Hosted Evo 2 model size: evo2-40b (default, most accurate) or evo2-7b (smaller, faster). One of `evo2-40b`, `evo2-7b` Default `evo2-40b` |
| `num_tokens` | integer | no | Number of nucleotides to generate Default `100` |
| `temperature` | number | no | Sampling temperature Default `0.7` |
| `top_k` | integer | no | Top-k sampling parameter Default `1` |
| `top_p` | number | no | Top-p (nucleus) sampling parameter |
| `enable_sampled_probs` | boolean | no | Return per-nucleotide probabilities Default `False` |
| `enable_logits` | boolean | no | Return raw logits Default `False` |

??? note "Upstream description"

    Generate DNA sequences using NVIDIA-hosted Evo 2 (Arc Institute genome foundation model, trained on 9 trillion nucleotides). Two hosted sizes are selectable via `model`: evo2-40b (default, most accurate) and evo2-7b (smaller, faster). Input: DNA sequence (ACTG), num_tokens to generate. Output: generated sequence with optional per-nucleotide probabilities. Fast synchronous response.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'NvidiaNIM_evo2', 'arguments': {'sequence': '...'}},
)
```

## `SpliceAI_get_max_delta` {#SpliceAI_get_max_delta}

Returns the maximum SpliceAI delta score across splice effects for a variant with a benign-to-high pathogenicity interpretation.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `variant`
- **Keywords:** splicing, delta score, splice variant, pathogenicity, acceptor donor
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant` | string | yes | Variant in chr-pos-ref-alt format (e.g., chr8-140300616-T-G) |
| `genome` | string | no | Genome build: 37 or 38 (default: 38) One of `37`, `38` Default `38` |

??? note "Upstream description"

    Get simplified maximum SpliceAI delta score and pathogenicity interpretation for a variant. Returns single max score across all splice effects (AG/AL/DG/DL) with interpretation: benign (&lt;0.2), low (0.2-0.5), moderate (0.5-0.8), high (≥0.8). Best for quick variant triage.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SpliceAI_get_max_delta', 'arguments': {'variant': '...'}},
)
```

## `SpliceAI_predict_pangolin` {#SpliceAI_predict_pangolin}

Predicts splice-altering effects of a variant per transcript with the Pangolin deep learning model.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `variant`
- **Keywords:** pangolin, splicing, splice gain loss, variant effect, deep learning
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant` | string | yes | Variant in chr-pos-ref-alt format (e.g., chr8-140300616-T-G) |
| `genome` | string | no | Genome build: 37 or 38 (default: 38) One of `37`, `38` Default `38` |
| `distance` | integer | no | Distance parameter for model (default: 1000) Default `1000` |
| `mask` | boolean | no | Use masked scores Default `False` |

??? note "Upstream description"

    Predict splice-altering effects using Pangolin model (alternative to SpliceAI). Pangolin uses a different architecture and may provide complementary predictions. Returns splice effect scores per transcript. Useful for validation or when SpliceAI gives borderline results.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SpliceAI_predict_pangolin',
     'arguments': {'variant': '...'}},
)
```

## `SpliceAI_predict_splice` {#SpliceAI_predict_splice}

Predicts splice-altering effects of a variant with SpliceAI, returning acceptor/donor gain and loss delta scores.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `variant`
- **Keywords:** splicing, acceptor gain, donor loss, delta score, variant interpretation, grch38
- **Read-only:** yes (declared by the server)

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `variant` | string | yes | Variant in chr-pos-ref-alt format (e.g., chr8-140300616-T-G) or colon-separated |
| `genome` | string | no | Genome build: 37 (GRCh37/hg19) or 38 (GRCh38/hg38). Default: 38 One of `37`, `38` Default `38` |
| `distance` | integer | no | Distance parameter for model (default: 50). Larger values capture more distal effects. Default `50` |
| `mask` | boolean | no | Use masked scores (recommended for variant interpretation). Raw scores better for alternative splicing analysis. Default `False` |

??? note "Upstream description"

    Predict splice-altering effects using SpliceAI deep learning model. Returns delta scores for acceptor gain (DS_AG), acceptor loss (DS_AL), donor gain (DS_DG), donor loss (DS_DL). Scores ≥0.2 suggest splice impact, ≥0.5 moderate, ≥0.8 high pathogenicity. Critical for intronic/exonic variant interpretation.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'tooluniverse',
    'execute_tool',
    {'tool_name': 'SpliceAI_predict_splice',
     'arguments': {'variant': '...'}},
)
```
