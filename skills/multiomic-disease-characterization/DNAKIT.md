# DNAKit integration

Canonical skill: `multiomic-disease-characterization`.

```bash
python -m pip install dnakit
dnakit skills doctor multiomic-disease-characterization
dnakit skills install multiomic-disease-characterization --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('clin_var_search_variants')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `ClinVar_search_variants` | `dnakit.databases.clin_var_search_variants` |
| `DGIdb_get_drug_gene_interactions` | `dnakit.databases.dgidb_get_drug_gene_interactions` |
| `ExpressionAtlas_search_differential` | `dnakit.databases.expression_atlas_search_differential` |
| `ExpressionAtlas_search_experiments` | `dnakit.databases.expression_atlas_search_experiments` |
| `GO_get_annotations_for_gene` | `dnakit.databases.go_get_annotations_for_gene` |
| `GO_search_terms` | `dnakit.databases.go_search_terms` |
| `GTEx_get_expression_summary` | `dnakit.databases.gtex_get_expression_summary` |
| `GWAS_search_associations_by_gene` | `dnakit.databases.gwas_search_associations_by_gene` |
| `HPA_get_cancer_prognostics_by_gene` | `dnakit.databases.hpa_get_cancer_prognostics_by_gene` |
| `HPA_get_protein_interactions_by_gene` | `dnakit.databases.hpa_get_protein_interactions_by_gene` |
| `HPA_get_rna_expression_by_source` | `dnakit.databases.hpa_get_rna_expression_by_source` |
| `HPA_get_rna_expression_in_specific_tissues` | `dnakit.databases.hpa_get_rna_expression_in_specific_tissues` |
| `HPA_get_subcellular_location` | `dnakit.databases.hpa_get_subcellular_location` |
| `HPA_search_genes_by_query` | `dnakit.databases.hpa_search_genes_by_query` |
| `MyGene_query_genes` | `dnakit.databases.my_gene_query_genes` |
| `OSL_get_efo_id_by_disease_name` | `dnakit.databases.osl_get_efo_id_by_disease_name` |
| `OpenTargets_get_associated_drugs_by_disease_efoId` | `dnakit.databases.open_targets_get_associated_drugs_by_disease_efo_id` |
| `OpenTargets_get_associated_drugs_by_target_ensemblID` | `dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id` |
| `OpenTargets_get_associated_phenotypes_by_disease_efoId` | `dnakit.databases.open_targets_get_associated_phenotypes_by_disease_efo_id` |
| `OpenTargets_get_associated_targets_by_disease_efoId` | `dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id` |
| `OpenTargets_get_disease_ancestors_parents_by_efoId` | `dnakit.databases.open_targets_get_disease_ancestors_parents_by_efo_id` |
| `OpenTargets_get_disease_descendants_children_by_efoId` | `dnakit.databases.open_targets_get_disease_descendants_children_by_efo_id` |
| `OpenTargets_get_disease_description_by_efoId` | `dnakit.databases.open_targets_get_disease_description_by_efo_id` |
| `OpenTargets_get_disease_id_description_by_name` | `dnakit.databases.open_targets_get_disease_id_description_by_name` |
| `OpenTargets_get_disease_synonyms_by_efoId` | `dnakit.databases.open_targets_get_disease_synonyms_by_efo_id` |
| `OpenTargets_get_disease_therapeutic_areas_by_efoId` | `dnakit.databases.open_targets_get_disease_therapeutic_areas_by_efo_id` |
| `OpenTargets_get_drug_mechanisms_of_action_by_chemblId` | `dnakit.databases.open_targets_get_drug_mechanisms_of_action_by_chembl_id` |
| `OpenTargets_get_evidence_by_datasource` | `dnakit.databases.open_targets_get_evidence_by_datasource` |
| `OpenTargets_get_similar_entities_by_disease_efoId` | `dnakit.databases.open_targets_get_similar_entities_by_disease_efo_id` |
| `OpenTargets_get_target_gene_ontology_by_ensemblID` | `dnakit.databases.open_targets_get_target_gene_ontology_by_ensembl_id` |
| `OpenTargets_get_target_interactions_by_ensemblID` | `dnakit.databases.open_targets_get_target_interactions_by_ensembl_id` |
| `OpenTargets_get_target_safety_profile_by_ensemblID` | `dnakit.databases.open_targets_get_target_safety_profile_by_ensembl_id` |
| `OpenTargets_get_target_tractability_by_ensemblID` | `dnakit.databases.open_targets_get_target_tractability_by_ensembl_id` |
| `OpenTargets_map_any_disease_id_to_all_other_ids` | `dnakit.databases.open_targets_map_any_disease_id_to_all_other_ids` |
| `OpenTargets_search_gwas_studies_by_disease` | `dnakit.databases.open_targets_search_gwas_studies_by_disease` |
| `PubMed_search_articles` | `dnakit.databases.pub_med_search_articles` |
| `QuickGO_annotations_by_gene` | `dnakit.databases.quick_go_annotations_by_gene` |
| `ReactomeAnalysis_pathway_enrichment` | `dnakit.databases.reactome_analysis_pathway_enrichment` |
| `Reactome_get_pathway` | `dnakit.databases.reactome_get_pathway` |
| `Reactome_get_pathway_reactions` | `dnakit.databases.reactome_get_pathway_reactions` |
| `Reactome_map_uniprot_to_pathways` | `dnakit.databases.reactome_map_uniprot_to_pathways` |
| `STRING_functional_enrichment` | `dnakit.databases.string_functional_enrichment` |
| `STRING_get_interaction_partners` | `dnakit.databases.string_get_interaction_partners` |
| `STRING_get_network` | `dnakit.databases.string_get_network` |
| `STRING_ppi_enrichment` | `dnakit.databases.string_ppi_enrichment` |
| `UniProt_get_function_by_accession` | `dnakit.databases.uni_prot_get_function_by_accession` |
| `WikiPathways_search` | `dnakit.databases.wiki_pathways_search` |
| `enrichr_gene_enrichment_analysis` | `dnakit.databases.enrichr_gene_enrichment_analysis` |
| `ensembl_lookup_gene` | `dnakit.databases.ensembl_lookup_gene` |
| `europepmc_disease_target_score` | `dnakit.databases.europepmc_disease_target_score` |
| `expression_atlas_disease_target_score` | `dnakit.databases.expression_atlas_disease_target_score` |
| `gnomad_get_gene_constraints` | `dnakit.databases.gnomad_get_gene_constraints` |
| `gwas_get_snps_for_gene` | `dnakit.databases.gwas_get_snps_for_gene` |
| `gwas_get_studies_for_trait` | `dnakit.databases.gwas_get_studies_for_trait` |
| `gwas_get_variants_for_trait` | `dnakit.databases.gwas_get_variants_for_trait` |
| `gwas_search_associations` | `dnakit.databases.gwas_search_associations` |
| `humanbase_ppi_analysis` | `dnakit.databases.humanbase_ppi_analysis` |
| `intact_get_interactions` | `dnakit.databases.intact_get_interactions` |
| `intact_search_interactions` | `dnakit.databases.intact_search_interactions` |
| `kegg_get_pathway_info` | `dnakit.databases.kegg_get_pathway_info` |
| `kegg_search_pathway` | `dnakit.databases.kegg_search_pathway` |
| `search_clinical_trials` | `dnakit.databases.search_clinical_trials` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
