# DNAKit integration

Canonical skill: `rare-disease-diagnosis`.

```bash
python -m pip install dnakit
dnakit skills doctor rare-disease-diagnosis
dnakit skills install rare-disease-diagnosis --agent codex
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
| `ArXiv_search_papers` | `dnakit.databases.ar_xiv_search_papers` |
| `BioRxiv_get_preprint` | `dnakit.databases.bio_rxiv_get_preprint` |
| `BioRxiv_list_recent_preprints` | `dnakit.databases.bio_rxiv_list_recent_preprints` |
| `CADD_get_variant_score` | `dnakit.databases.cadd_get_variant_score` |
| `CELLxGENE_download_h5ad` | `dnakit.databases.cellx_gene_download_h5ad` |
| `CELLxGENE_get_cell_metadata` | `dnakit.databases.cellx_gene_get_cell_metadata` |
| `CELLxGENE_get_embeddings` | `dnakit.databases.cellx_gene_get_embeddings` |
| `CELLxGENE_get_expression_data` | `dnakit.databases.cellx_gene_get_expression_data` |
| `ChIPAtlas_enrichment_analysis` | `dnakit.databases.ch_ipatlas_enrichment_analysis` |
| `ChIPAtlas_get_experiments` | `dnakit.databases.ch_ipatlas_get_experiments` |
| `ChIPAtlas_get_peak_data` | `dnakit.databases.ch_ipatlas_get_peak_data` |
| `ChIPAtlas_search_datasets` | `dnakit.databases.ch_ipatlas_search_datasets` |
| `ClinGen_get_variant_classifications` | `dnakit.databases.clin_gen_get_variant_classifications` |
| `ClinGen_search_actionability` | `dnakit.databases.clin_gen_search_actionability` |
| `ClinGen_search_dosage_sensitivity` | `dnakit.databases.clin_gen_search_dosage_sensitivity` |
| `ClinGen_search_gene_validity` | `dnakit.databases.clin_gen_search_gene_validity` |
| `ClinVar_get_clinical_significance` | `dnakit.databases.clin_var_get_clinical_significance` |
| `ClinVar_get_variant_details` | `dnakit.databases.clin_var_get_variant_details` |
| `ClinVar_search_variants` | `dnakit.databases.clin_var_search_variants` |
| `Crossref_search_works` | `dnakit.databases.crossref_search_works` |
| `DisGeNET_get_disease_genes` | `dnakit.databases.dis_ge_net_get_disease_genes` |
| `DisGeNET_get_gda` | `dnakit.databases.dis_ge_net_get_gda` |
| `DisGeNET_get_vda` | `dnakit.databases.dis_ge_net_get_vda` |
| `DisGeNET_search_disease` | `dnakit.databases.dis_ge_net_search_disease` |
| `DisGeNET_search_gene` | `dnakit.databases.dis_ge_net_search_gene` |
| `ENCODE_get_biosample` | `dnakit.databases.encode_get_biosample` |
| `ENCODE_get_experiment` | `dnakit.databases.encode_get_experiment` |
| `ENCODE_list_files` | `dnakit.databases.encode_list_files` |
| `ENCODE_search_experiments` | `dnakit.databases.encode_search_experiments` |
| `ESM_describe_sae_feature` | `dnakit.databases.esm_describe_sae_feature` |
| `ESM_explain_variant_mechanism` | `dnakit.databases.esm_explain_variant_mechanism` |
| `ESM_get_region_sae_features` | `dnakit.databases.esm_get_region_sae_features` |
| `ESM_score_variant_sae_batch` | `dnakit.databases.esm_score_variant_sae_batch` |
| `ESM_score_variant_sae_disruption` | `dnakit.databases.esm_score_variant_sae_disruption` |
| `EVE_get_variant_score` | `dnakit.databases.eve_get_variant_score` |
| `EuropePMC_search_articles` | `dnakit.databases.europe_pmc_search_articles` |
| `FAVOR_annotate_variant` | `dnakit.databases.favor_annotate_variant` |
| `GTEx_get_expression_summary` | `dnakit.databases.gtex_get_expression_summary` |
| `GTEx_get_median_gene_expression` | `dnakit.databases.gtex_get_median_gene_expression` |
| `HPA_get_rna_expression_by_source` | `dnakit.databases.hpa_get_rna_expression_by_source` |
| `HPA_search_genes_by_query` | `dnakit.databases.hpa_search_genes_by_query` |
| `HPO_get_diseases_by_phenotype` | `dnakit.databases.hpo_get_diseases_by_phenotype` |
| `HPO_get_genes_by_phenotype` | `dnakit.databases.hpo_get_genes_by_phenotype` |
| `HPO_get_term` | `dnakit.databases.hpo_get_term` |
| `HPO_search_terms` | `dnakit.databases.hpo_search_terms` |
| `InterPro_get_protein_domains` | `dnakit.databases.inter_pro_get_protein_domains` |
| `KEGG_get_gene_pathways` | `dnakit.databases.kegg_get_gene_pathways` |
| `MARRVEL_get_gene` | `dnakit.databases.marrvel_get_gene` |
| `MARRVEL_get_omim_phenotypes` | `dnakit.databases.marrvel_get_omim_phenotypes` |
| `MyGene_get_gene_annotation` | `dnakit.databases.my_gene_get_gene_annotation` |
| `MyGene_query_genes` | `dnakit.databases.my_gene_query_genes` |
| `NvidiaNIM_alphafold2` | `dnakit.databases.nvidia_nim_alphafold2` |
| `NvidiaNIM_esmfold` | `dnakit.databases.nvidia_nim_esmfold` |
| `OMIM_get_clinical_synopsis` | `dnakit.databases.omim_get_clinical_synopsis` |
| `OMIM_get_entry` | `dnakit.databases.omim_get_entry` |
| `OMIM_get_gene_map` | `dnakit.databases.omim_get_gene_map` |
| `OMIM_search` | `dnakit.databases.omim_search` |
| `OpenTargets_get_associated_drugs_by_target_ensemblID` | `dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id` |
| `OpenTargets_get_associated_targets_by_disease_efoId` | `dnakit.databases.open_targets_get_associated_targets_by_disease_efo_id` |
| `OpenTargets_get_disease_description_by_efoId` | `dnakit.databases.open_targets_get_disease_description_by_efo_id` |
| `OpenTargets_get_diseases_phenotypes_by_target_ensembl` | `dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl` |
| `Orphanet_get_classification` | `dnakit.databases.orphanet_get_classification` |
| `Orphanet_get_disease` | `dnakit.databases.orphanet_get_disease` |
| `Orphanet_get_genes` | `dnakit.databases.orphanet_get_genes` |
| `Orphanet_search_by_name` | `dnakit.databases.orphanet_search_by_name` |
| `Orphanet_search_diseases` | `dnakit.databases.orphanet_search_diseases` |
| `Pfam_get_protein_annotations` | `dnakit.databases.pfam_get_protein_annotations` |
| `PubMed_get_article` | `dnakit.databases.pub_med_get_article` |
| `PubMed_get_cited_by` | `dnakit.databases.pub_med_get_cited_by` |
| `PubMed_get_related` | `dnakit.databases.pub_med_get_related` |
| `PubMed_search_articles` | `dnakit.databases.pub_med_search_articles` |
| `ReactomeContent_search` | `dnakit.databases.reactome_content_search` |
| `Reactome_get_pathway` | `dnakit.databases.reactome_get_pathway` |
| `STRING_get_interaction_partners` | `dnakit.databases.string_get_interaction_partners` |
| `SemanticScholar_search_papers` | `dnakit.databases.semantic_scholar_search_papers` |
| `SpliceAI_get_max_delta` | `dnakit.databases.splice_ai_get_max_delta` |
| `SpliceAI_predict_pangolin` | `dnakit.databases.splice_ai_predict_pangolin` |
| `SpliceAI_predict_splice` | `dnakit.databases.splice_ai_predict_splice` |
| `UniProt_get_entry_by_accession` | `dnakit.databases.uni_prot_get_entry_by_accession` |
| `UniProt_get_features_by_accession` | `dnakit.databases.uni_prot_get_features_by_accession` |
| `alphafold_get_prediction` | `dnakit.databases.alphafold_get_prediction` |
| `ensembl_lookup_gene` | `dnakit.databases.ensembl_lookup_gene` |
| `get_joint_associated_diseases_by_HPO_ID_list` | `dnakit.databases.get_joint_associated_diseases_by_hpo_id_list` |
| `gnomad_get_gene_constraints` | `dnakit.databases.gnomad_get_gene_constraints` |
| `gnomad_get_variant` | `dnakit.databases.gnomad_get_variant` |
| `intact_get_complex_details` | `dnakit.databases.intact_get_complex_details` |
| `intact_get_interaction_network` | `dnakit.databases.intact_get_interaction_network` |
| `intact_search_interactions` | `dnakit.databases.intact_search_interactions` |
| `kegg_find_genes` | `dnakit.databases.kegg_find_genes` |
| `kegg_get_gene_info` | `dnakit.databases.kegg_get_gene_info` |
| `kegg_get_pathway_info` | `dnakit.databases.kegg_get_pathway_info` |
| `kegg_search_pathway` | `dnakit.databases.kegg_search_pathway` |
| `openalex_get_author` | `dnakit.databases.openalex_get_author` |
| `openalex_literature_search` | `dnakit.databases.openalex_literature_search` |
| `openalex_search_works` | `dnakit.databases.openalex_search_works` |
| `reactome_disease_target_score` | `dnakit.databases.reactome_disease_target_score` |
| `web_search` | `dnakit.databases.web_search` |

## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts rare-disease-diagnosis
dnakit skills run rare-disease-diagnosis --script scripts/clinical_patterns.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `scripts/clinical_patterns.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.

## Regulatory result types

`ch_ipatlas_get_peak_data(gene=...)` resolves the GRCh38 gene locus and obtains an assembled track link. A link response has `status="requires_data"`; it contains no observed peaks. Download and verify the selected BED, then pass `peak_file=...` to compute interval overlaps. `flank_bp` and `limit` control the queried interval and returned rows. This is not an antigen-name search.

The inherited `ch_ipatlas_enrichment_analysis` adapter supplies web-submission instructions, not computed enrichment scores; DNAKit marks this `requires_submission`. Complete that external analysis and inspect actual statistics before interpreting TF binding. Peak overlap alone does not establish motif disruption. See [ChIP-Atlas's API](https://chip-atlas.org/agents).
