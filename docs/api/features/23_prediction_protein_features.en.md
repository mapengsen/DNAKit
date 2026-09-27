# Protein Features and Structure

See the [FAQ](../../faq/tutorial_integrations.md#property-prediction) for dependencies, checkpoints and runtime configuration.

<!-- BEGIN GENERATED API: model-proteins -->

<span id="model-proteins"></span>

## 1) Label biological categories of protein SAE features {#dnakit.databases.esm_describe_sae_feature}

- **Purpose:** Label biological categories of protein SAE features.
- **Method:** Use `ESM` to run the selected trained model on the supplied input.
- **API：** `dnakit.databases.esm_describe_sae_feature(*, feature_id, sae_model=None, model=None, n_proteins=None, top_residues_per_protein=None, use_cache=None, verbose=True)`
- **Input:** Provide `feature_id`. Configure: `ESM_API_KEY`.
- **Output:** Returns a mapping with `feature_id`, `sae_model`, `category`, `confidence` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.esm_describe_sae_feature(
    feature_id=6404,
    n_proteins=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "feature_id": "<integer>",
  "sae_model": "<string>",
  "category": "<string>",
  "confidence": "<number>",
  "n_proteins_with_activation": "<integer>",
  "n_proteins_analyzed": "<integer>",
  "category_vote_counts": "<object>",
  "supporting_evidence": [
    {
      "protein": "<string>",
      "position_1based": "<integer>",
      "activation": "<number>",
      "uniprot_types": "<array>",
      "uniprot_categories": "<array>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_protein_features.md#dnakit.databases.esm_describe_sae_feature-parameters)

## 2) Aggregate protein SAE features over a residue interval {#dnakit.databases.esm_get_region_sae_features}

- **Purpose:** Aggregate protein SAE features over a residue interval.
- **Method:** Use `ESM` to run the selected trained model on the supplied input.
- **API：** `dnakit.databases.esm_get_region_sae_features(*, sequence, start_position, end_position, top_k_features=None, model=None, sae_model=None, verbose=True)`
- **Input:** Provide `sequence`, `start_position`, `end_position`. Configure: `ESM_API_KEY`.
- **Output:** Returns a mapping with `sequence_length`, `region`, `region_length`, `n_features_active_in_region` and related fields.
- **Example code:**

```python
import dnakit

example_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEK'

result = dnakit.databases.esm_get_region_sae_features(
    sequence=example_sequence,
    start_position=10,
    end_position=17,
    top_k_features=10,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "sequence_length": "<integer>",
  "region": [
    "<value>"
  ],
  "region_length": "<integer>",
  "n_features_active_in_region": "<integer>",
  "top_features": [
    "<value>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_protein_features.md#dnakit.databases.esm_get_region_sae_features-parameters)

## 3) Extract per-residue protein SAE activations {#dnakit.databases.esm_get_sae_features}

- **Purpose:** Extract per-residue protein SAE activations.
- **Method:** Use `ESM` to run the selected trained model on the supplied input.
- **API：** `dnakit.databases.esm_get_sae_features(*, sequence, model=None, sae_model=None, position=None, window=None, top_k_per_residue=None, verbose=True)`
- **Input:** Provide `sequence`. Configure: `ESM_API_KEY`.
- **Output:** Returns a mapping with `sequence_length`, `model`, `sae_model`, `residues_returned` and related fields.
- **Example code:**

```python
import dnakit

example_sequence = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDR'

result = dnakit.databases.esm_get_sae_features(
    sequence=example_sequence,
    position=175,
    window=8,
    top_k_per_residue=16,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "sequence_length": "<integer>",
  "model": "<string>",
  "sae_model": "<string>",
  "residues_returned": "<integer>",
  "position": "<integer>",
  "window": "<integer>",
  "top_k_per_residue": "<integer>",
  "activations": [
    {
      "residue_idx_1based": "<integer>",
      "active_features": "<array>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_protein_features.md#dnakit.databases.esm_get_sae_features-parameters)

## 4) Score protein sequence log-probabilities with ESMC {#dnakit.databases.esm_score_sequence}

- **Purpose:** Score protein sequence log-probabilities with ESMC.
- **Method:** Use `ESM` to run the selected trained model on the supplied input.
- **API：** `dnakit.databases.esm_score_sequence(*, sequence, model=None, verbose=True)`
- **Input:** Provide `sequence`. Configure: `ESM_API_KEY`.
- **Output:** Returns a mapping with `status`, `model`, `sequence`, `sequence_length` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.esm_score_sequence(
    sequence='MKTAYIAKQRQISFVKSHFSRQ',
    model='esmc-300m-2024-12',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "status": "<string>",
  "model": "<string>",
  "sequence": "<string>",
  "sequence_length": "<integer>",
  "mean_log_probability": "<number>",
  "per_residue_log_probabilities": [
    "<number>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_protein_features.md#dnakit.databases.esm_score_sequence-parameters)

## 5) Predict protein structures with AlphaFold2 {#dnakit.databases.nvidia_nim_alphafold2}

- **Purpose:** Predict protein structures with AlphaFold2.
- **Method:** Use `NvidiaNIM` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.nvidia_nim_alphafold2(*, sequence, algorithm=None, databases=None, e_value=None, iterations=None, relax_prediction=None, skip_template_search=None, verbose=True)`
- **Input:** Provide `sequence`. Configure: `NVIDIA_API_KEY`.
- **Output:** Returns predicted PDB coordinates and service-provided confidence and model metrics.
- **Example code:**

```python
import dnakit

result = dnakit.databases.nvidia_nim_alphafold2(
    sequence='MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSH',
    verbose=False,
)
print(result)
```

- **Example output:** Returns predicted PDB coordinates and service-provided confidence and model metrics.

[Full parameters and usage notes](../../faq/api/23_prediction_protein_features.md#dnakit.databases.nvidia_nim_alphafold2-parameters)

## 6) Predict protein structures from sequence with ESMFold {#dnakit.databases.nvidia_nim_esmfold}

- **Purpose:** Predict protein structures from sequence with ESMFold.
- **Method:** Use `NvidiaNIM` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.nvidia_nim_esmfold(*, sequence, verbose=True)`
- **Input:** Provide `sequence`. Configure: `NVIDIA_API_KEY`.
- **Output:** Returns a string.
- **Example code:**

```python
import dnakit

result = dnakit.databases.nvidia_nim_esmfold(
    sequence='ARNDCQEGHILKMFPSTWYV',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
"<string>"
```

[Full parameters and usage notes](../../faq/api/23_prediction_protein_features.md#dnakit.databases.nvidia_nim_esmfold-parameters)

<!-- END GENERATED API: model-proteins -->
