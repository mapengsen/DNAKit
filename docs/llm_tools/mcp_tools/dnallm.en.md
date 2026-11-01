# DNALLM MCP Server (`dnallm`)

DNAKit's MCP tool index keeps **11** DNA/gene tools of this server (16 listed, read on 2026-10-16). They run once you install the server; DNAKit does not bundle them. See [calling a tool that was found](../../api/features/111_mcp_tool_index.md#4-calling-a-tool-that-was-found).

| Item | Value |
| --- | --- |
| Source | [zhangtaolab/DNALLM @ `f7a7fcb`](https://github.com/zhangtaolab/DNALLM/tree/f7a7fcb914f94a2d82c2af96be0e648eab4a719f) |
| License | MIT |
| Command | `dnallm-mcp-server` |
| Transports | stdio |
| Installation | Install the PyPI package dnallm (Python &gt;=3.11, installs the dnallm-mcp-server console script), plus a torch build suited to the hardware; because the wheel does not ship dnallm/mcp/configs, also keep a source checkout and set workdir to it so the default --config dnallm/mcp/configs/mcp_server_config.yaml resolves. |
| Effects | Loads the DNA language models named in the YAML configuration (downloading weights from Hugging Face or ModelScope on first use), runs local CPU/GPU inference, writes temporary VCF files for variant scoring, and may save predictions under the output_dir set in the model inference configs; no sequence data is sent to a DNALLM-operated service. |

| Tool | Category | Summary |
| --- | --- | --- |
| [`_dna_batch_predict`](#_dna_batch_predict) | `prediction_model` | Predicts properties of multiple DNA sequences in one batch using a single loaded DNA language model. |
| [`_dna_interpret`](#_dna_interpret) | `prediction_model` | Computes per-token attribution scores for a DNA sequence prediction using Captum methods such as integrated gradients, DeepLift or occlusion. |
| [`_dna_multi_model_predict`](#_dna_multi_model_predict) | `prediction_model` | Predicts one DNA sequence with multiple loaded models in parallel for model comparison or ensembling. |
| [`_dna_mutagenesis`](#_dna_mutagenesis) | `prediction_model` | Performs in silico mutagenesis (substitutions, insertions, deletions) on DNA sequences and scores the effect on model predictions. |
| [`_dna_sequence_predict`](#_dna_sequence_predict) | `prediction_model` | Predicts properties of a single DNA sequence using a specified loaded DNA model. |
| [`_dna_stream_batch_predict`](#_dna_stream_batch_predict) | `prediction_model` | Predicts multiple DNA sequences with one model while streaming progress updates. |
| [`_dna_stream_multi_model_predict`](#_dna_stream_multi_model_predict) | `prediction_model` | Predicts one DNA sequence with multiple models while streaming progress updates. |
| [`_dna_stream_predict`](#_dna_stream_predict) | `prediction_model` | Predicts one DNA sequence with a model while streaming progress updates via server-sent events. |
| [`_hotspots`](#_hotspots) | `prediction_model` | Identifies mutational hotspot windows in a genomic region (up to 2 kb from a FASTA) by model in silico mutagenesis and sliding-window aggregation. |
| [`_ism_scan`](#_ism_scan) | `prediction_model` | Runs a bounded in silico mutagenesis scan (up to 2000 bases, 100 reported positions) on DNA sequences with a loaded model. |
| [`_zero_shot_score`](#_zero_shot_score) | `prediction_model` | Scores zero-shot variant effects with a DNA model for inline chrom/pos/ref/alt variants or a VCF, using a reference FASTA and optional ClinVar filter. |

## `_dna_batch_predict` {#_dna_batch_predict}

Predicts properties of multiple DNA sequences in one batch using a single loaded DNA language model.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequences`, `model_name`
- **Keywords:** dna language model, batch inference, sequence classification, genomic model

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequences` | array of string | yes |  |
| `model_name` | string | yes |  |

??? note "Upstream description"

    Predict multiple DNA sequences using a specific model.

    This tool performs batch prediction on multiple DNA sequences using
    a single model. It's optimized for processing multiple sequences
    efficiently by leveraging model batching capabilities.

    Args:
    sequences (list\[str\]): List of DNA sequences to predict.
    Each sequence should contain only valid nucleotide
    characters (A, T, G, C). Case insensitive.
    model_name (str): Name of the model to use for all predictions.
    Must be one of the loaded models.

    Returns:
    dict\[str, Any\]: Batch prediction results in MCP format:
    - On success: Contains 'content', 'model_name',
    'sequence_count'
    - On error: Contains 'error', 'isError' fields

    Example:
    ```python
    result = await server._dna_batch_predict(
    sequences=\["ATCGATCG", "GCTAGCTA", "TTAACCGG"\],
    model_name=(
    "zhangtaolab/plant-dnabert-BPE-open_chromatin"
    )
    )
    ```

    Note:
    Batch processing is generally more efficient than individual
    predictions for multiple sequences, especially with GPU
    acceleration.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'dnallm',
    '_dna_batch_predict',
    {'sequences': [], 'model_name': '...'},
)
```

## `_dna_interpret` {#_dna_interpret}

Computes per-token attribution scores for a DNA sequence prediction using Captum methods such as integrated gradients, DeepLift or occlusion.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequence`, `model_name`
- **Keywords:** interpretability, attribution, captum, integrated gradients, saliency, dna language model

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes |  |
| `model_name` | string | yes |  |
| `method` | string | no | Default `lig` |
| `target_class` | integer \| null | no | Default `None` |
| `max_length` | integer \| null | no | Default `None` |

??? note "Upstream description"

    Interpret model predictions using attribution methods.

    This tool provides model interpretability by computing attribution
    scores for each token in a DNA sequence using various Captum methods.

    Args:
    sequence (str): DNA sequence to interpret.
    model_name (str): Name of the model to use.
    method (str): Attribution method. One of: "lig",
    "deeplift", "occlusion", "feature_ablation",
    "layer_conductance", "gradient_shap", "noise_tunnel",
    "integrated_gradients". Defaults to "lig".
    target_class (int \| None): Target class index for attribution.
    If None, auto-selects the class with maximum probability.
    max_length (int \| None): Maximum token length for tokenizer.

    Returns:
    dict\[str, Any\]: Interpretation results in MCP format:
    - On success: Contains 'content', 'attributions',
    'tokens', 'method', 'target_class', 'model_name',
    'sequence'
    - On error: Contains 'error', 'isError' fields

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'dnallm',
    '_dna_interpret',
    {'sequence': '...', 'model_name': '...'},
)
```

## `_dna_multi_model_predict` {#_dna_multi_model_predict}

Predicts one DNA sequence with multiple loaded models in parallel for model comparison or ensembling.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequence`
- **Keywords:** ensemble, model comparison, dna language model, consensus

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes |  |
| `model_names` | array \| null | no | Default `None` |

??? note "Upstream description"

    Predict DNA sequence using multiple models in parallel.

    This tool performs prediction on a single DNA sequence using multiple
    models simultaneously, enabling model comparison and ensemble analysis.
    It's particularly useful for understanding prediction consensus across
    different model architectures.

    Args:
    sequence (str): DNA sequence to predict, containing only valid
    nucleotide characters (A, T, G, C). Case insensitive.
    model_names (list\[str\] \| None, optional): List of model names
    to use. If None, uses all currently loaded models.
    Defaults to None.

    Returns:
    dict\[str, Any\]: Multi-model prediction results in MCP format:
    - On success: Contains 'content', 'model_count', 'sequence'
    - On error: Contains 'error', 'isError' fields

    Example:
    ```python
    \# Use specific models
    result = await server._dna_multi_model_predict(
    sequence="ATCGATCG",
    model_names=\[
    "dnabert-2", "nucleotide-transformer"
    \]
    )

    \# Use all loaded models
    result = await server._dna_multi_model_predict(
    sequence="ATCGATCG"
    )
    ```

    Note:
    Models are processed in parallel for better performance.
    Individual model failures don't stop the entire operation.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'dnallm',
    '_dna_multi_model_predict',
    {'sequence': '...'},
)
```

## `_dna_mutagenesis` {#_dna_mutagenesis}

Performs in silico mutagenesis (substitutions, insertions, deletions) on DNA sequences and scores the effect on model predictions.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `model_name`
- **Keywords:** in silico mutagenesis, ism, mutation effect, substitution, indel, dna language model

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `model_name` | string | yes |  |
| `sequence` | string \| null | no | Default `None` |
| `sequences` | array \| null | no | Default `None` |
| `mutation_type` | string | no | Default `single_base_substitution` |
| `positions` | array \| null | no | Default `None` |

??? note "Upstream description"

    Perform in silico mutagenesis on DNA sequences.

    This tool evaluates the impact of sequence mutations on model
    predictions, supporting single base substitutions, multi-base
    substitutions, deletions, insertions, and exhaustive combinations.

    Args:
    sequence (str \| None): Single DNA sequence to mutate. If provided,
    it is processed as a one-element list internally.
    sequences (list\[str\] \| None): List of DNA sequences to mutate.
    Either sequence or sequences must be provided.
    mutation_type (str): Type of mutation to perform. One of:
    "single_base_substitution", "multi_base_substitution",
    "deletion", "insertion", "combo". Defaults to
    "single_base_substitution".
    positions (list\[int\] \| None): 0-based positions to mutate.
    Required and must be non-empty.
    model_name (str): Name of the model to use for prediction.

    Returns:
    dict\[str, Any\]: Mutagenesis results in MCP format:
    - On success: Contains 'content', 'original_prediction',
    'mutated_prediction', 'delta', 'affected_positions',
    'mutation_type', 'model_name'
    - On error: Contains 'error', 'isError' fields

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'dnallm',
    '_dna_mutagenesis',
    {'model_name': '...'},
)
```

## `_dna_sequence_predict` {#_dna_sequence_predict}

Predicts properties of a single DNA sequence using a specified loaded DNA model.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequence`, `model_name`
- **Keywords:** dna language model, inference, sequence classification, genomic model

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes |  |
| `model_name` | string | yes |  |

??? note "Upstream description"

    Predict DNA sequence using a specific model.

    This tool performs single DNA sequence prediction using a specified
    pre-loaded model. It's the most basic prediction operation and serves
    as the foundation for more complex prediction tasks.

    Args:
    sequence (str): DNA sequence to predict, containing only valid
    nucleotide characters (A, T, G, C). Case insensitive.
    model_name (str): Name of the model to use for prediction.
    Must be one of the loaded models.

    Returns:
    dict\[str, Any\]: Prediction results in MCP format:
    - On success: Contains 'content', 'model_name', 'sequence'
    - On error: Contains 'error', 'isError' fields

    Example:
    ```python
    result = await server._dna_sequence_predict(
    sequence="ATCGATCG",
    model_name="dnabert-2"
    )
    ```

    Note:
    This method handles all exceptions internally and returns
    error information in the response rather than raising exceptions.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'dnallm',
    '_dna_sequence_predict',
    {'sequence': '...', 'model_name': '...'},
)
```

## `_dna_stream_batch_predict` {#_dna_stream_batch_predict}

Predicts multiple DNA sequences with one model while streaming progress updates.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequences`, `model_name`
- **Keywords:** streaming, batch inference, progress, dna language model

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequences` | array of string | yes |  |
| `model_name` | string | yes |  |
| `stream_progress` | boolean | no | Default `True` |
| `context` | ? \| null | no | Default `None` |

??? note "Upstream description"

    Stream batch DNA sequence prediction with real-time progress
    updates.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'dnallm',
    '_dna_stream_batch_predict',
    {'sequences': [], 'model_name': '...'},
)
```

## `_dna_stream_multi_model_predict` {#_dna_stream_multi_model_predict}

Predicts one DNA sequence with multiple models while streaming progress updates.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequence`
- **Keywords:** streaming, ensemble, progress, dna language model

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes |  |
| `model_names` | array \| null | no | Default `None` |
| `stream_progress` | boolean | no | Default `True` |
| `context` | ? \| null | no | Default `None` |

??? note "Upstream description"

    Stream multi-model DNA sequence prediction with real-time
    progress updates.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'dnallm',
    '_dna_stream_multi_model_predict',
    {'sequence': '...'},
)
```

## `_dna_stream_predict` {#_dna_stream_predict}

Predicts one DNA sequence with a model while streaming progress updates via server-sent events.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `sequence`, `model_name`
- **Keywords:** streaming, sse, inference, dna language model

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `sequence` | string | yes |  |
| `model_name` | string | yes |  |
| `stream_progress` | boolean | no | Default `True` |
| `context` | ? \| null | no | Default `None` |

??? note "Upstream description"

    Stream DNA sequence prediction with real-time progress updates.

    This tool provides real-time streaming prediction with progress updates
    via Server-Sent Events (SSE). It's designed for interactive
    applications
    where users need immediate feedback on prediction progress.

    The streaming capability is particularly useful for:
    - Long-running predictions on large sequences
    - Interactive web applications requiring real-time feedback
    - Progress monitoring for batch operations

    Args:
    sequence (str): DNA sequence to predict, containing only valid
    nucleotide characters (A, T, G, C). Case insensitive.
    model_name (str): Name of the model to use for prediction.
    Must be one of the loaded models.
    stream_progress (bool, optional): Whether to send progress updates
    via SSE. Defaults to True.
    context (Any \| None, optional): MCP Context object for progress
    reporting. Required for progress updates. Defaults to None.

    Returns:
    dict\[str, Any\]: Streaming prediction results in MCP format:
    - On success: Contains 'content', 'model_name',
    'sequence', 'streamed'
    - On error: Contains 'error', 'isError' fields

    Example:
    ```python
    \# With progress streaming
    result = await server._dna_stream_predict(
    sequence="ATCGATCG",
    model_name="dnabert-2",
    stream_progress=True,
    context=mcp_context
    )
    ```

    Note:
    Progress updates are sent at key stages: initialization (0%),
    model loading (25%), processing (75%), and completion (100%).
    Uses chunk-based timeout that resets after each progress update.

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'dnallm',
    '_dna_stream_predict',
    {'sequence': '...', 'model_name': '...'},
)
```

## `_hotspots` {#_hotspots}

Identifies mutational hotspot windows in a genomic region (up to 2 kb from a FASTA) by model in silico mutagenesis and sliding-window aggregation.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `model_name`, `coordinates`, `fasta_path`
- **Keywords:** mutational hotspot, in silico mutagenesis, sliding window, genomic region, fasta

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `model_name` | string | yes |  |
| `coordinates` | object | yes |  |
| `fasta_path` | string | yes |  |
| `strategy` | string | no | Default `maxabs` |
| `window_size` | integer | no | Default `10` |
| `percentile_threshold` | number | no | Default `90.0` |

??? note "Upstream description"

    Identify mutational hotspot windows for a genomic region.

    Computes hotspots from the model itself (in silico mutagenesis over
    the requested reference slice followed by sliding-window
    aggregation) — never from a precomputed window file. The reference
    sequence is read from a per-call server-side ``fasta_path``.

    Args:
    model_name (str): Name of the loaded model to scan with.
    coordinates (dict\[str, Any\]): Region to scan, as
    ``\{"chrom": str, "start": int, "end": int\}`` with 0-based
    half-open coordinates (Python slicing convention). The
    region length is capped at 2000 bases per call.
    fasta_path (str): Server-side reference FASTA path
    (.fasta/.fa/.fa.gz/.fna) containing the chromosome. This is
    an operator-trust-boundary file read, not a client upload.
    strategy (str): Per-base score aggregation, one of "maxabs",
    "min", "max", "mean". Defaults to "maxabs".
    window_size (int): Sliding-window size (bases) for hotspot
    detection. Defaults to 10.
    percentile_threshold (float): Rolling-window percentile above
    which a window is a hotspot (0 &lt; p &lt;= 100). Defaults to 90.

    Returns:
    dict\[str, Any\]: Hotspot scan results in MCP format:
    - On success: Contains 'content', 'hotspots' (0-based
    half-open \[start, end) pairs relative to the region),
    'hotspots_genomic' (absolute coordinates),
    'window_count', 'coordinates', 'sequence_length',
    'strategy', 'window_size', 'percentile_threshold',
    'model_name'
    - On error: Contains 'error', 'isError' fields

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'dnallm',
    '_hotspots',
    {'model_name': '...', 'coordinates': {}, 'fasta_path': '...'},
)
```

## `_ism_scan` {#_ism_scan}

Runs a bounded in silico mutagenesis scan (up to 2000 bases, 100 reported positions) on DNA sequences with a loaded model.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `model_name`
- **Keywords:** in silico mutagenesis, ism, saturation mutagenesis, mutation effect, dna language model

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `model_name` | string | yes |  |
| `sequence` | string \| null | no | Default `None` |
| `sequences` | array \| null | no | Default `None` |
| `mutation_type` | string | no | Default `single_base_substitution` |
| `positions` | array \| null | no | Default `None` |

??? note "Upstream description"

    Run a bounded in silico mutagenesis (ISM) scan on DNA sequences.

    Mirrors the ``dna_mutagenesis`` engine surface with tool-boundary
    input caps so a scan cannot silently expand past the tool timeout:
    sequences are capped at 2000 bases and the requested position list
    at 100 entries per call (single-base-substitution ISM costs 3
    forward passes per base). The engine itself scans every position of
    every sequence — the bound on COMPUTE is the sequence cap — while
    the REPORTED mutated entries are restricted to the requested
    positions. Reduce the inputs or raise ``tool_timeout_seconds`` in
    the server config when a larger scan is genuinely needed.

    Args:
    model_name (str): Name of the loaded model to scan with.
    sequence (str \| None): Single DNA sequence to scan. Mutually
    exclusive with ``sequences`` (passing both is a matchable
    error; the single sequence is otherwise processed as a
    one-element list internally).
    sequences (list\[str\] \| None): Non-empty list of non-empty DNA
    sequence strings to scan. Mutually exclusive with
    ``sequence``; an empty list is rejected like missing input.
    mutation_type (str): One of "single_base_substitution",
    "multi_base_substitution", "deletion", "insertion", "combo".
    positions (list\[int\] \| None): 0-based positions of interest;
    validated against every sequence's length and echoed in the
    response. Must be non-empty. The reported mutated entries
    (and the delta averages) are restricted to these positions;
    the engine still scans every position of every sequence, so
    the bound on compute is the sequence cap, not this list.

    Returns:
    dict\[str, Any\]: ISM scan results in MCP format:
    - On success: Contains 'content', 'original_prediction',
    'mutated_prediction', 'delta', 'affected_positions',
    'mutation_type', 'model_name'
    - On error: Contains 'error', 'isError' fields

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'dnallm',
    '_ism_scan',
    {'model_name': '...'},
)
```

## `_zero_shot_score` {#_zero_shot_score}

Scores zero-shot variant effects with a DNA model for inline chrom/pos/ref/alt variants or a VCF, using a reference FASTA and optional ClinVar filter.

- **Category:** `prediction_model` (DNA sequence models and variant-effect prediction models)
- **Required arguments:** `model_name`, `fasta_path`
- **Keywords:** zero-shot, variant effect prediction, vcf, clinvar, log-likelihood ratio, dna language model

| Argument | Type | Required | Description |
| --- | --- | --- | --- |
| `model_name` | string | yes |  |
| `fasta_path` | string | yes |  |
| `variants` | array \| null | no | Default `None` |
| `vcf_path` | string \| null | no | Default `None` |
| `paradigm` | string | no | Default `mlm` |
| `context_window` | integer | no | Default `200` |
| `clnsig_filter` | object \| null | no | Default `None` |

??? note "Upstream description"

    Score zero-shot variant effects (dual-mode, D-04).

    Accepts EITHER an inline ``variants`` list of \{chrom, pos, ref, alt\}
    objects (pos is the 1-based VCF coordinate; materialized to a
    server-side temp VCF with a fixed sanitized name) OR a server-side
    ``vcf_path`` (.vcf/.vcf.gz, size-capped). Both modes route through
    the SAME ``dnallm.inference.vep.evaluate_vcf`` kernel, so skip
    accounting (skip_counts, skipped, skip_fraction) and the convention
    block are identical and surfaced verbatim in the response.

    Inline variants carry no ClinVar annotation, so they are admitted
    under a pass-through convention (labels "not_analyzed"=1, metrics
    None by construction — see the convention block in the response).
    Pass ``clnsig_filter`` (fields: variant_type, positive_labels,
    negative_labels, star_floor — mapping to ``vep.ClinVarFilter``) to
    override the D-17 defaults, e.g. for non-ClinVar ``vcf_path`` input.
    In inline-variants mode the pass-through sentinel is the mandatory
    base: ``variant_type`` is always ``inline_variant`` (an explicit
    conflicting value is a matchable error naming the required value)
    and the caller's other fields override the sentinel, not the D-17
    defaults — any other ``variant_type`` would exclude every inline
    row as ``non_snv_clnvc``.

    Variant count is capped at 500 per call; sequences are capped at
    2000 bases for the ISM tools. ``fasta_path``/``vcf_path`` are
    operator-trust-boundary server-side file reads, never uploads.

    Args:
    model_name (str): Name of the loaded model to score with.
    fasta_path (str): Server-side reference FASTA path
    (.fasta/.fa/.fa.gz/.fna) for window building.
    variants (list\[dict\[str, Any\]\] \| None): Inline variants, each
    \{"chrom": str (whitespace-free chromosome name), "pos": int
    (1-based), "ref": ACGT str, "alt": ACGT str\}. Mutually
    exclusive with vcf_path.
    vcf_path (str \| None): Server-side ClinVar-style VCF path.
    Mutually exclusive with variants.
    paradigm (str): "mlm" (log-odds, default) or "clm"
    (delta-log-likelihood).
    context_window (int): Reference bases kept on each side of a
    variant. Defaults to 200.
    clnsig_filter (dict\[str, Any\] \| None): Optional ClinVar
    convention override (see above); ``None`` keeps the D-17
    defaults for vcf_path mode and the pass-through convention
    for inline mode. In inline mode the sentinel
    ``variant_type`` is mandatory (see above).

    Returns:
    dict\[str, Any\]: Zero-shot scoring results in MCP format:
    - On success: Contains 'content', 'records' (per-variant
    scores/skips), 'skip_counts', 'evaluated', 'skipped',
    'skip_fraction', 'metrics', 'convention' (all verbatim
    from the kernel's VepResult), 'input_mode', 'paradigm',
    'model_name'
    - On error: Contains 'error', 'isError' fields

Call it through DNAKit once the server is installed:

```python
from dnakit.integrations import call_mcp_tool

result = call_mcp_tool(
    'dnallm',
    '_zero_shot_score',
    {'model_name': '...', 'fasta_path': '...'},
)
```
