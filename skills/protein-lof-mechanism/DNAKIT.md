# DNAKit integration

Canonical skill: `protein-lof-mechanism`.

```bash
python -m pip install dnakit
dnakit skills doctor protein-lof-mechanism
dnakit skills install protein-lof-mechanism --agent codex
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
| `DynaMut2_predict_stability` | `dnakit.databases.dyna_mut2_predict_stability` |
| `ESM_describe_sae_feature` | `dnakit.databases.esm_describe_sae_feature` |
| `ESM_explain_variant_mechanism` | `dnakit.databases.esm_explain_variant_mechanism` |
| `ESM_score_sequence` | `dnakit.databases.esm_score_sequence` |
| `ESM_score_variant_sae_disruption` | `dnakit.databases.esm_score_variant_sae_disruption` |
| `UniProt_get_entry_by_accession` | `dnakit.databases.uni_prot_get_entry_by_accession` |
| `UniProt_get_sequence_by_accession` | `dnakit.databases.uni_prot_get_sequence_by_accession` |
| `UniProt_search` | `dnakit.databases.uni_prot_search` |
| `alphafold_get_prediction` | `dnakit.databases.alphafold_get_prediction` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
