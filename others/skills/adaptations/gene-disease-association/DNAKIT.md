# DNAKit integration

Canonical skill: `gene-disease-association`.

```bash
python -m pip install dnakit
dnakit skills doctor gene-disease-association
dnakit skills install gene-disease-association --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('clin_var_get_variant_details')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `ClinVar_get_variant_details` | `dnakit.databases.clin_var_get_variant_details` |
| `ClinVar_search_variants` | `dnakit.databases.clin_var_search_variants` |
| `DisGeNET_get_disease_genes` | `dnakit.databases.dis_ge_net_get_disease_genes` |
| `DisGeNET_get_gda` | `dnakit.databases.dis_ge_net_get_gda` |
| `DisGeNET_get_vda` | `dnakit.databases.dis_ge_net_get_vda` |
| `DisGeNET_search_disease` | `dnakit.databases.dis_ge_net_search_disease` |
| `DisGeNET_search_gene` | `dnakit.databases.dis_ge_net_search_gene` |
| `GenCC_get_classifications` | `dnakit.databases.gen_cc_get_classifications` |
| `GenCC_search_disease` | `dnakit.databases.gen_cc_search_disease` |
| `GenCC_search_gene` | `dnakit.databases.gen_cc_search_gene` |
| `Harmonizome_get_gene` | `dnakit.databases.harmonizome_get_gene` |
| `MonarchV3_get_associations` | `dnakit.databases.monarch_v3_get_associations` |
| `MonarchV3_get_entity` | `dnakit.databases.monarch_v3_get_entity` |
| `MonarchV3_get_histopheno` | `dnakit.databases.monarch_v3_get_histopheno` |
| `MonarchV3_get_mappings` | `dnakit.databases.monarch_v3_get_mappings` |
| `MonarchV3_search` | `dnakit.databases.monarch_v3_search` |
| `Monarch_search` | `dnakit.databases.monarch_search` |
| `MyGene_query_genes` | `dnakit.databases.my_gene_query_genes` |
| `OMIM_get_clinical_synopsis` | `dnakit.databases.omim_get_clinical_synopsis` |
| `OMIM_get_entry` | `dnakit.databases.omim_get_entry` |
| `OMIM_get_gene_map` | `dnakit.databases.omim_get_gene_map` |
| `OMIM_search` | `dnakit.databases.omim_search` |
| `OpenTargets_get_diseases_phenotypes_by_target_ensembl` | `dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl` |
| `OpenTargets_multi_entity_search_by_query_string` | `dnakit.databases.open_targets_multi_entity_search_by_query_string` |
| `OpenTargets_target_disease_evidence` | `dnakit.databases.open_targets_target_disease_evidence` |
| `Orphanet_get_gene_diseases` | `dnakit.databases.orphanet_get_gene_diseases` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
