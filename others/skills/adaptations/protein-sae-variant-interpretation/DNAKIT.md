# DNAKit integration

Canonical skill: `protein-sae-variant-interpretation`.

```bash
python -m pip install dnakit
dnakit skills doctor protein-sae-variant-interpretation
dnakit skills install protein-sae-variant-interpretation --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('alpha_missense_get_variant_score')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `AlphaMissense_get_variant_score` | `dnakit.databases.alpha_missense_get_variant_score` |
| `ClinVar_get_variant_details` | `dnakit.databases.clin_var_get_variant_details` |
| `ClinVar_search_variants` | `dnakit.databases.clin_var_search_variants` |
| `ESM_describe_sae_feature` | `dnakit.databases.esm_describe_sae_feature` |
| `ESM_explain_variant_mechanism` | `dnakit.databases.esm_explain_variant_mechanism` |
| `ESM_get_sae_features` | `dnakit.databases.esm_get_sae_features` |
| `ESM_score_variant_sae_batch` | `dnakit.databases.esm_score_variant_sae_batch` |
| `ESM_score_variant_sae_disruption` | `dnakit.databases.esm_score_variant_sae_disruption` |
| `UniProt_get_disease_variants_by_accession` | `dnakit.databases.uni_prot_get_disease_variants_by_accession` |
| `UniProt_get_function_by_accession` | `dnakit.databases.uni_prot_get_function_by_accession` |
| `UniProt_get_sequence_by_accession` | `dnakit.databases.uni_prot_get_sequence_by_accession` |
| `UniProt_search` | `dnakit.databases.uni_prot_search` |
| `alphafold_get_prediction` | `dnakit.databases.alphafold_get_prediction` |
| `gnomad_get_variant` | `dnakit.databases.gnomad_get_variant` |
| `gnomad_search_variants` | `dnakit.databases.gnomad_search_variants` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
