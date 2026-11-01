# Enformer predictions and metrics

## Purpose

Explicit DNA operations for integration pack P39/P40, with optional upstream backends and JSON results.

## Method

Each function invokes the named upstream operation. Local paths are explicit; new outputs must not already exist. UCSC/blacklist retrieval functions use the network only when called. NUPACK and BLAST+ are not present in the validation environment; those scientific pipelines remain unvalidated. Enformer checkpoints must be local and match the supplied architecture and adapter configuration.

## API

`dnakit.predictions.enformer_p40_functions`

| Function | Operation |
| --- | --- |
| `enformer_poisson_loss` | Compute upstream poisson_loss, including its epsilon and mean reduction. |
| `enformer_pearson_correlation` | Compute upstream per-sequence Pearson coefficient for batch x positions x tracks. |
| `enformer_channel_correlation` | Accumulate upstream MeanPearsonCorrCoefPerChannel over successive track batches. |
| `enformer_predict_tracks` | Predict the human or mouse tracks using a complete local Enformer checkpoint. |
| `enformer_predict_custom_tracks` | Predict user-defined tracks with a trained HeadAdapterWrapper checkpoint. |
| `enformer_predict_context_tracks` | Predict context-specific tracks using trained hypernetwork adapter weights. |
| `enformer_predict_context_attention` | Predict tracks with trained cross-attention context-token adapter weights. |

## Input

DNA strings, JSON mappings/lists, local genomic files or explicit local model checkpoints. Interval coordinates and physical-model parameters follow each function docstring; the NUPACK default material is DNA. All functions accept keyword-only `verbose`.

## Output

A dictionary with completion status, backend provenance and operation-specific data. Writers return output paths. Figures return image paths and underlying data. Optional backends raise `BackendUnavailableError` when missing. P40 manifests distinguish unvalidated pipelines from tested interfaces.

**Example code:**

```python
from dnakit.predictions.enformer_p40_functions import enformer_poisson_loss
print(enformer_poisson_loss([1.0], [1.0], verbose=False)['values'])
```

**Example output:**

```text
1.0
```
