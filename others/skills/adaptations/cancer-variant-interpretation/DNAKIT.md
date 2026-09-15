# DNAKit integration

Canonical skill: `cancer-variant-interpretation`.

```bash
python -m pip install dnakit
dnakit skills doctor cancer-variant-interpretation
dnakit skills install cancer-variant-interpretation --agent codex
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
| `ChEMBL_search_drugs` | `dnakit.databases.ch_embl_search_drugs` |
| `ESM_explain_variant_mechanism` | `dnakit.databases.esm_explain_variant_mechanism` |
| `FDA_get_boxed_warning_info_by_drug_name` | `dnakit.databases.fda_get_boxed_warning_info_by_drug_name` |
| `FDA_get_clinical_studies_info_by_drug_name` | `dnakit.databases.fda_get_clinical_studies_info_by_drug_name` |
| `FDA_get_indications_by_drug_name` | `dnakit.databases.fda_get_indications_by_drug_name` |
| `FDA_get_mechanism_of_action_by_drug_name` | `dnakit.databases.fda_get_mechanism_of_action_by_drug_name` |
| `GTEx_get_median_gene_expression` | `dnakit.databases.gtex_get_median_gene_expression` |
| `MyGene_query_genes` | `dnakit.databases.my_gene_query_genes` |
| `OpenTargets_get_associated_drugs_by_disease_efoId` | `dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id` |
| `OpenTargets_get_associated_drugs_by_target_ensemblID` | `dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id` |
| `OpenTargets_get_disease_id_description_by_name` | `dnakit.databases.open_targets_get_disease_id_description_by_name` |
| `OpenTargets_get_drug_adverse_events_by_chemblId` | `dnakit.databases.open_targets_get_drug_adverse_events_by_chembl_id` |
| `OpenTargets_get_drug_chembId_by_generic_name` | `dnakit.databases.open_targets_get_drug_chemb_id_by_generic_name` |
| `OpenTargets_get_drug_indications_by_chemblId` | `dnakit.databases.open_targets_get_drug_indications_by_chembl_id` |
| `OpenTargets_get_drug_mechanisms_of_action_by_chemblId` | `dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id` |
| `OpenTargets_get_publications_by_target_ensemblID` | `dnakit.databases.open_targets_get_publications_by_target_ensembl_id` |
| `OpenTargets_get_target_id_description_by_name` | `dnakit.databases.open_targets_get_target_id_description_by_name` |
| `OpenTargets_target_disease_evidence` | `dnakit.databases.open_targets_target_disease_evidence` |
| `PubMed_search_articles` | `dnakit.databases.pub_med_search_articles` |
| `Reactome_map_uniprot_to_pathways` | `dnakit.databases.reactome_map_uniprot_to_pathways` |
| `UniProt_get_disease_variants_by_accession` | `dnakit.databases.uni_prot_get_disease_variants_by_accession` |
| `UniProt_get_function_by_accession` | `dnakit.databases.uni_prot_get_function_by_accession` |
| `UniProt_search` | `dnakit.databases.uni_prot_search` |
| `cBioPortal_get_cancer_studies` | `dnakit.databases.c_bio_portal_get_cancer_studies` |
| `cBioPortal_get_gene_info` | `dnakit.databases.c_bio_portal_get_gene_info` |
| `cBioPortal_get_molecular_profiles` | `dnakit.databases.c_bio_portal_get_molecular_profiles` |
| `cBioPortal_get_mutations` | `dnakit.databases.c_bio_portal_get_mutations` |
| `cBioPortal_get_samples` | `dnakit.databases.c_bio_portal_get_samples` |
| `civic_get_molecular_profile` | `dnakit.databases.civic_get_molecular_profile` |
| `civic_get_variant` | `dnakit.databases.civic_get_variant` |
| `civic_get_variants_by_gene` | `dnakit.databases.civic_get_variants_by_gene` |
| `civic_search_assertions` | `dnakit.databases.civic_search_assertions` |
| `civic_search_evidence_items` | `dnakit.databases.civic_search_evidence_items` |
| `civic_search_genes` | `dnakit.databases.civic_search_genes` |
| `civic_search_therapies` | `dnakit.databases.civic_search_therapies` |
| `drugbank_get_drug_basic_info_by_drug_name_or_id` | `dnakit.databases.drugbank_get_drug_basic_info_by_drug_name_or_id` |
| `drugbank_get_pharmacology_by_drug_name_or_drugbank_id` | `dnakit.databases.drugbank_get_pharmacology_by_drug_name_or_drugbank_id` |
| `drugbank_get_targets_by_drug_name_or_drugbank_id` | `dnakit.databases.drugbank_get_targets_by_drug_name_or_drugbank_id` |
| `ensembl_lookup_gene` | `dnakit.databases.ensembl_lookup_gene` |
| `search_clinical_trials` | `dnakit.databases.search_clinical_trials` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
