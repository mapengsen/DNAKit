# DNAKit integration

Canonical skill: `variant-predictor-dms-validation`.

```bash
python -m pip install dnakit
dnakit skills doctor variant-predictor-dms-validation
dnakit skills install variant-predictor-dms-validation --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('alpha_missense_get_protein_scores')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `AlphaMissense_get_protein_scores` | `dnakit.databases.alpha_missense_get_protein_scores` |
| `AlphaMissense_get_residue_scores` | `dnakit.databases.alpha_missense_get_residue_scores` |
| `AlphaMissense_get_variant_score` | `dnakit.databases.alpha_missense_get_variant_score` |
| `ESM2_score_missense_variant` | `dnakit.databases.esm2_score_missense_variant` |
| `ESM_get_sae_features` | `dnakit.databases.esm_get_sae_features` |
| `ESM_score_sequence` | `dnakit.databases.esm_score_sequence` |
| `ESM_score_variant_sae_batch` | `dnakit.databases.esm_score_variant_sae_batch` |
| `ESM_score_variant_sae_disruption` | `dnakit.databases.esm_score_variant_sae_disruption` |
| `MaveDB_get_effect_matrix` | `dnakit.databases.mave_db_get_effect_matrix` |
| `Structure_annotate_per_residue` | `dnakit.databases.structure_annotate_per_residue` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
