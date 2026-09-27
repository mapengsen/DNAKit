# Variant Effects and Protein Stability

See the [FAQ](../../faq/tutorial_integrations.md#property-prediction) for dependencies, checkpoints and runtime configuration.

## 1) Long-context zero-shot variant effects {#15-long-context-zero-shot-variant-effects}

- **Purpose:** Calculate REF and ALT sequence likelihoods and an ALT−REF effect score for a single SNV.
- **Prediction method:** Use the official Evo 2 7B base model to calculate mean sequence log-likelihood over the complete REF and ALT contexts, then compute `alternate − reference`. This is a zero-shot ranking score, not a clinically calibrated pathogenicity probability.
- **API:** `dnakit.predictions.predict_variant_effects(inputs[required], config=PropertyPredictionConfig(model="evo2", task="variant_effect")[required], backend[optional])`
- **Input:** `VariantContext`; REF and ALT must have equal lengths and differ at exactly one A/C/G/T position.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

```python
from dnakit.predictions import (
    PropertyPredictionConfig,
    VariantContext,
    predict_variant_effects,
)

variant = VariantContext(
    "snv-1",
    reference_sequence="A" * 4096 + "C" + "G" * 4095,
    alternate_sequence="A" * 4096 + "T" + "G" * 4095,
)
result = predict_variant_effects(
    [variant],
    config=PropertyPredictionConfig(
        model="evo2",
        task="variant_effect",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.values.shape)
```

- **Example output:**

```text
('reference_likelihood', 'alternate_likelihood', 'alternate_minus_reference')
(3,)
```

## 2) Allele-conditional variant effects {#17-allele-conditional-variant-effects}

- **Purpose:** Compare next-token support for REF and ALT bases from the sequence upstream of a variant.
- **Prediction method:** Use the official GENERator v2 eukaryote 1.2B weights with up to 8,192 bp of upstream context. Remove leading `N`, left-truncate to a multiple of 6, and calculate `log(p_ref/(p_alt+1e-10))`. This is a zero-shot ranking score, not a clinically calibrated probability.
- **API:** `dnakit.predictions.predict_variant_effects(inputs[required], config=PropertyPredictionConfig(model="generator", task="variant_effect", allow_remote_code=True)[required], backend[optional])`
- **Input:** `VariantContext`; REF and ALT must have equal lengths and exactly one SNV. Only the reference context before the variant is used.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

```python
from dnakit.predictions import (
    PropertyPredictionConfig,
    VariantContext,
    predict_variant_effects,
)

variant = VariantContext(
    "snv-1",
    reference_sequence="ACGT" * 1024 + "C" + "G" * 255,
    alternate_sequence="ACGT" * 1024 + "T" + "G" * 255,
)
result = predict_variant_effects(
    [variant],
    config=PropertyPredictionConfig(
        model="generator",
        task="variant_effect",
        allow_remote_code=True,
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.values.shape)
```

- **Example output:**

```text
('reference_probability', 'alternate_probability', 'log_ratio_score')
(3,)
```

## 3) Protein stability {#23-protein-stability}

- **Purpose:** Predict a stability regression value from a protein sequence under the official ProtStab task definition.
- **Prediction method:** Run the official LucaOneTasks `ProtStab` regression checkpoint through `predict_v1.py`. DNAKit does not fit a regressor.
- **API:** `dnakit.predictions.predict_sequence_properties(inputs[required], config=PropertyPredictionConfig(model="lucaone", task="protein_stability", model_source_path=...)[required], backend[optional])`
- **Input:** `BiologicalSequence` with `sequence_type="protein"`.
- **Output:** Returns `PropertyPredictionResult`; each entry in `records` has prediction `output_names`, numeric `values`, and related metadata.
- **Example code:**

```python
from dnakit.predictions import (
    BiologicalSequence,
    PropertyPredictionConfig,
    predict_sequence_properties,
)

result = predict_sequence_properties(
    [
        BiologicalSequence(
            "protein-1",
            "MAMAPRTEINSTRING" * 10,
            "protein",
        )
    ],
    config=PropertyPredictionConfig(
        model="lucaone",
        task="protein_stability",
        model_source_path="/opt/LucaOneTasks",
    ),
)
output = result.records[0].output
print(output.output_names)
print(output.values.shape)
```

- **Example output (value depends on input and checkpoint):**

```text
('prediction',)
(1,)
```

<!-- BEGIN GENERATED API: model-variants -->

<span id="model-variants"></span>

## 4) Score regulatory variant effects with AlphaGenome {#dnakit.databases.alpha_genome_score_variant}

- **Purpose:** Score regulatory variant effects with AlphaGenome.
- **Method:** Use `AlphaGenome` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.alpha_genome_score_variant(*, chromosome, position, reference_bases, alternate_bases, output_type=None, organism=None, sequence_length=None, top_n=None, verbose=True)`
- **Input:** Provide `chromosome`, `position`, `reference_bases`, `alternate_bases`. Configure: `ALPHA_GENOME_API_KEY`.
- **Output:** Returns a mapping with `variant`, `output_type`, `scores` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.alpha_genome_score_variant(
    chromosome='chr22',
    position=36201698,
    reference_bases='A',
    alternate_bases='C',
    output_type='RNA_SEQ',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "variant": "<string>",
  "output_type": "<string>",
  "scores": [
    "<object>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.alpha_genome_score_variant-parameters)

## 5) Predict mutation effects on protein stability with DynaMut2 {#dnakit.databases.dyna_mut2_predict_stability}

- **Purpose:** Predict mutation effects on protein stability with DynaMut2.
- **Method:** Use `DynaMut2` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.dyna_mut2_predict_stability(*, operation=None, pdb_id, chain, mutation, verbose=True)`
- **Input:** Provide `pdb_id`, `chain`, `mutation`.
- **Output:** Returns a mapping with `pdb_id`, `chain`, `mutation`, `wild_type` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.dyna_mut2_predict_stability(
    operation='predict_stability',
    pdb_id='4HHB',
    chain='A',
    mutation='V1A',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "pdb_id": "<string>",
  "chain": "<string>",
  "mutation": "<string>",
  "wild_type": "<string>",
  "mutant": "<string>",
  "position": "<string>",
  "ddg_prediction": "<number>",
  "effect": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.dyna_mut2_predict_stability-parameters)

## 6) Score missense variants with ESM-2 {#dnakit.databases.esm2_score_missense_variant}

- **Purpose:** Score missense variants with ESM-2.
- **Method:** Use `ESM2VariantEffect` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.esm2_score_missense_variant(*, sequence, position, mutant, wild_type=None, model_id=None, wait_for_model=None, verbose=True)`
- **Input:** Provide `sequence`, `position`, `mutant`.
- **Output:** Returns a mapping with `variant`, `wild_type`, `mutant`, `p_wild_type` and related fields.
- **Example code:**

```python
import dnakit

example_sequence = 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'

result = dnakit.databases.esm2_score_missense_variant(
    sequence=example_sequence,
    position=7,
    wild_type='E',
    mutant='V',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "variant": "<string>",
  "wild_type": "<string>",
  "mutant": "<string>",
  "p_wild_type": "<number>",
  "p_mutant": "<number>",
  "log_likelihood_ratio": "<number>",
  "direction": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.esm2_score_missense_variant-parameters)

## 7) Summarize protein feature disruption caused by a variant {#dnakit.databases.esm_explain_variant_mechanism}

- **Purpose:** Summarize protein feature disruption caused by a variant.
- **Method:** Use `ESM` to run the selected trained model on the supplied input.
- **API：** `dnakit.databases.esm_explain_variant_mechanism(*, sequence, position, ref_aa, alt_aa, window=None, top_k_features=None, include_descriptions=None, model=None, sae_model=None, verbose=True)`
- **Input:** Provide `sequence`, `position`, `ref_aa`, `alt_aa`. Configure: `ESM_API_KEY`.
- **Output:** Returns a mapping with `variant`, `position`, `window`, `mechanism_summary` and related fields.
- **Example code:**

```python
import dnakit

example_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEK'

result = dnakit.databases.esm_explain_variant_mechanism(
    sequence=example_sequence,
    position=12,
    ref_aa='G',
    alt_aa='D',
    top_k_features=3,
    include_descriptions=False,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "variant": "<string>",
  "position": "<integer>",
  "window": "<integer>",
  "mechanism_summary": "<string>",
  "lost_feature_categories": "<object>",
  "gained_feature_categories": "<object>",
  "top_features_lost": [
    "<value>"
  ],
  "top_features_gained": [
    "<value>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.esm_explain_variant_mechanism-parameters)

## 8) Batch-score missense variants using SAE disruption {#dnakit.databases.esm_score_variant_sae_batch}

- **Purpose:** Batch-score missense variants using SAE disruption.
- **Method:** Use `ESM` to run the selected trained model on the supplied input.
- **API：** `dnakit.databases.esm_score_variant_sae_batch(*, sequence, variants, window=None, top_k_features=None, model=None, sae_model=None, verbose=True)`
- **Input:** Provide `sequence`, `variants`. Configure: `ESM_API_KEY`.
- **Output:** Returns a mapping with `sequence_length`, `n_variants`, `n_succeeded`, `window` and related fields.
- **Example code:**

```python
import dnakit

example_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEK'

example_variants = [{'position': 12, 'ref_aa': 'G', 'alt_aa': 'D'},
 {'position': 12, 'ref_aa': 'G', 'alt_aa': 'V'},
 {'position': 11, 'ref_aa': 'A', 'alt_aa': 'T'}]

result = dnakit.databases.esm_score_variant_sae_batch(
    sequence=example_sequence,
    variants=example_variants,
    top_k_features=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "sequence_length": "<integer>",
  "n_variants": "<integer>",
  "n_succeeded": "<integer>",
  "window": "<integer>",
  "results": [
    "<value>"
  ]
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.esm_score_variant_sae_batch-parameters)

## 9) Compare reference and mutant protein SAE features {#dnakit.databases.esm_score_variant_sae_disruption}

- **Purpose:** Compare reference and mutant protein SAE features.
- **Method:** Use `ESM` to run the selected trained model on the supplied input.
- **API：** `dnakit.databases.esm_score_variant_sae_disruption(*, sequence, position, ref_aa, alt_aa, window=None, top_k_features=None, model=None, sae_model=None, verbose=True)`
- **Input:** Provide `sequence`, `position`, `ref_aa`, `alt_aa`. Configure: `ESM_API_KEY`.
- **Output:** Returns a mapping with `variant`, `position`, `window`, `n_unique_features_touched` and related fields.
- **Example code:**

```python
import dnakit

example_sequence = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDR'

result = dnakit.databases.esm_score_variant_sae_disruption(
    sequence=example_sequence,
    position=175,
    ref_aa='R',
    alt_aa='H',
    window=8,
    top_k_features=5,
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "variant": "<string>",
  "position": "<integer>",
  "window": "<integer>",
  "n_unique_features_touched": "<integer>",
  "top_features_lost": [
    {
      "feature_id": "<integer>",
      "delta": "<number>",
      "ref_activation_sum": "<number>",
      "mut_activation_sum": "<number>"
    }
  ],
  "top_features_gained": [
    {
      "feature_id": "<integer>",
      "delta": "<number>",
      "ref_activation_sum": "<number>",
      "mut_activation_sum": "<number>"
    }
  ]
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.esm_score_variant_sae_disruption-parameters)

## 10) Score DNA variant effects with Evo 2 {#dnakit.databases.evo2_score_variant}

- **Purpose:** Score DNA variant effects with Evo 2.
- **Method:** Use `Evo2VariantEffect` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.evo2_score_variant(*, ref_sequence=None, alt_sequence=None, sequence=None, position=None, reference=None, alternate=None, model=None, verbose=True)`
- **Input:** Choose filters as shown in the example; omitted options use the operation defaults. Configure: `NVIDIA_API_KEY`.
- **Output:** Returns a mapping with `delta_loglik`, `ref_loglik`, `alt_loglik`, `direction` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.evo2_score_variant(
    sequence='ACGTACGTACGTACGTACGTGGGCCCAAATTTACGTACGTACGT',
    position=21,
    reference='G',
    alternate='T',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "delta_loglik": "<number>",
  "ref_loglik": "<number>",
  "alt_loglik": "<number>",
  "direction": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.evo2_score_variant-parameters)

<!-- END GENERATED API: model-variants -->

<!-- BEGIN GENERATED API: variant-splicing -->

<span id="variant-splicing"></span>

## 11) Retrieve a variant's maximum SpliceAI delta score {#dnakit.databases.splice_ai_get_max_delta}

- **Purpose:** Retrieve a variant's maximum SpliceAI delta score.
- **Method:** Use `SpliceAI` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.splice_ai_get_max_delta(*, variant, genome=None, verbose=True)`
- **Input:** Provide `variant`.
- **Output:** Returns a mapping with `data`, `source` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.splice_ai_get_max_delta(
    variant='chr8-140300616-T-G',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "variant": "<string>",
    "genome": "<string>",
    "max_delta_score": "<number>",
    "interpretation": "<string>",
    "pathogenicity_threshold": "<string>"
  },
  "source": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.splice_ai_get_max_delta-parameters)

## 12) Predict splice effects with Pangolin {#dnakit.databases.splice_ai_predict_pangolin}

- **Purpose:** Predict splice effects with Pangolin.
- **Method:** Use `SpliceAI` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.splice_ai_predict_pangolin(*, variant, genome=None, distance=None, mask=None, verbose=True)`
- **Input:** Provide `variant`.
- **Output:** Returns a mapping with `data`, `source` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.splice_ai_predict_pangolin(
    variant='chr8-140300616-T-G',
    genome='38',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "variant": "<string>",
    "genome": "<string>",
    "scores": [
      "<value>"
    ]
  },
  "source": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.splice_ai_predict_pangolin-parameters)

## 13) Predict splice-site changes with SpliceAI {#dnakit.databases.splice_ai_predict_splice}

- **Purpose:** Predict splice-site changes with SpliceAI.
- **Method:** Use `SpliceAI` to query or submit the supplied input and parse the response.
- **API：** `dnakit.databases.splice_ai_predict_splice(*, variant, genome=None, distance=None, mask=None, verbose=True)`
- **Input:** Provide `variant`.
- **Output:** Returns a mapping with `data`, `source` and related fields.
- **Example code:**

```python
import dnakit

result = dnakit.databases.splice_ai_predict_splice(
    variant='chr8-140300616-T-G',
    genome='38',
    verbose=False,
)
print(result)
```

- **Example output:** The response has the following shape; records depend on the input and source.

Response-shape excerpt; `<…>` denotes type placeholders, not a captured query result.

```json
{
  "data": {
    "variant": "<string>",
    "genome": "<string>",
    "scores": [
      "<value>"
    ],
    "max_delta_score": "<number>",
    "interpretation": "<string>"
  },
  "source": "<string>"
}
```

[Full parameters and usage notes](../../faq/api/23_prediction_variant_effects.md#dnakit.databases.splice_ai_predict_splice-parameters)

<!-- END GENERATED API: variant-splicing -->
