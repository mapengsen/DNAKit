# DNAKit integration

Canonical skill: `gwas-drug-discovery`.

```bash
python -m pip install dnakit
dnakit skills doctor gwas-drug-discovery
dnakit skills install gwas-drug-discovery --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('ch_embl_get_drug_mechanisms')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `ChEMBL_get_drug_mechanisms` | `dnakit.databases.ch_embl_get_drug_mechanisms` |
| `ChEMBL_get_molecule` | `dnakit.databases.ch_embl_get_molecule` |
| `ChEMBL_get_target_activities` | `dnakit.databases.ch_embl_get_target_activities` |
| `ChEMBL_search_drugs` | `dnakit.databases.ch_embl_search_drugs` |
| `ChEMBL_search_targets` | `dnakit.databases.ch_embl_search_targets` |
| `ClinicalTrials_search_studies` | `dnakit.databases.clinical_trials_search_studies` |
| `DGIdb_get_drug_gene_interactions` | `dnakit.databases.dgidb_get_drug_gene_interactions` |
| `EuropePMC_search_articles` | `dnakit.databases.europe_pmc_search_articles` |
| `FDA_get_active_ingredient_info_by_drug_name` | `dnakit.databases.fda_get_active_ingredient_info_by_drug_name` |
| `FDA_get_adverse_reactions_by_drug_name` | `dnakit.databases.fda_get_adverse_reactions_by_drug_name` |
| `OpenTargets_get_associated_drugs_by_disease_efoId` | `dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id` |
| `OpenTargets_get_associated_targets_by_disease_efoId` | `dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id` |
| `OpenTargets_get_disease_id_description_by_name` | `dnakit.databases.open_targets_get_disease_id_description_by_name` |
| `OpenTargets_get_drug_chembId_by_generic_name` | `dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name` |
| `OpenTargets_get_drug_mechanisms_of_action_by_chemblId` | `dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id` |
| `OpenTargets_get_drug_warnings_by_chemblId` | `dnakit.databases.open_targets_get_drug_warnings_by_chembl_id` |
| `OpenTargets_get_target_classes_by_ensemblID` | `dnakit.databases.open_targets_get_target_classes_by_ensembl_id` |
| `OpenTargets_get_target_genomic_location_by_ensemblID` | `dnakit.databases.open_targets_get_target_genomic_location_by_ensembl_id` |
| `OpenTargets_get_target_safety_profile_by_ensemblID` | `dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id` |
| `OpenTargets_get_target_tractability_by_ensemblID` | `dnakit.databases.open_targets_get_target_tractability_by_ensembl_id` |
| `OpenTargets_get_variant_credible_sets` | `dnakit.databases.open_targets_get_variant_credible_sets` |
| `OpenTargets_get_variant_info` | `dnakit.databases.open_targets_get_variant_info` |
| `OpenTargets_search_gwas_studies_by_disease` | `dnakit.databases.open_targets_search_gwas_studies_by_disease` |
| `PubMed_search_articles` | `dnakit.databases.pub_med_search_articles` |
| `execute_tool` | `dnakit.databases.execute_tool` |
| `gwas_get_associations_for_snp` | `dnakit.databases.gwas_get_associations_for_snp` |
| `gwas_get_associations_for_trait` | `dnakit.databases.gwas_get_associations_for_trait` |
| `gwas_search_associations` | `dnakit.databases.gwas_search_associations` |
| `gwas_search_snps` | `dnakit.databases.gwas_search_snps` |
| `gwas_search_studies` | `dnakit.databases.gwas_search_studies` |

## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts gwas-drug-discovery
dnakit skills run gwas-drug-discovery --script python_implementation.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `python_implementation.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
