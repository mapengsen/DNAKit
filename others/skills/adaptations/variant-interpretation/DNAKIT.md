# DNAKit integration

Canonical skill: `variant-interpretation`.

```bash
python -m pip install dnakit
dnakit skills doctor variant-interpretation
dnakit skills install variant-interpretation --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('alpha_genome_score_variant')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `AlphaGenome_score_variant` | `dnakit.databases.alpha_genome_score_variant` |
| `AlphaMissense_get_residue_scores` | `dnakit.databases.alpha_missense_get_residue_scores` |
| `AlphaMissense_get_variant_score` | `dnakit.databases.alpha_missense_get_variant_score` |
| `BioRxiv_get_preprint` | `dnakit.databases.bio_rxiv_get_preprint` |
| `BioRxiv_list_recent_preprints` | `dnakit.databases.bio_rxiv_list_recent_preprints` |
| `CADD_get_position_scores` | `dnakit.databases.cadd_get_position_scores` |
| `CADD_get_range_scores` | `dnakit.databases.cadd_get_range_scores` |
| `CADD_get_variant_score` | `dnakit.databases.cadd_get_variant_score` |
| `CELLxGENE_get_cell_metadata` | `dnakit.databases.cellx_gene_get_cell_metadata` |
| `CELLxGENE_get_expression_data` | `dnakit.databases.cellx_gene_get_expression_data` |
| `COSMIC_get_mutations_by_gene` | `dnakit.databases.cosmic_get_mutations_by_gene` |
| `COSMIC_search_mutations` | `dnakit.databases.cosmic_search_mutations` |
| `ChIPAtlas_enrichment_analysis` | `dnakit.databases.ch_ipatlas_enrichment_analysis` |
| `ChIPAtlas_get_peak_data` | `dnakit.databases.ch_ipatlas_get_peak_data` |
| `ChIPAtlas_search_datasets` | `dnakit.databases.ch_ipatlas_search_datasets` |
| `ClinGen_dosage_by_gene` | `dnakit.databases.clin_gen_dosage_by_gene` |
| `ClinGen_get_actionability_adult` | `dnakit.databases.clin_gen_get_actionability_adult` |
| `ClinGen_get_actionability_pediatric` | `dnakit.databases.clin_gen_get_actionability_pediatric` |
| `ClinGen_get_dosage_sensitivity` | `dnakit.databases.clin_gen_get_dosage_sensitivity` |
| `ClinGen_get_gene_validity` | `dnakit.databases.clin_gen_get_gene_validity` |
| `ClinGen_search_actionability` | `dnakit.databases.clin_gen_search_actionability` |
| `ClinGen_search_dosage_sensitivity` | `dnakit.databases.clin_gen_search_dosage_sensitivity` |
| `ClinGen_search_gene_validity` | `dnakit.databases.clin_gen_search_gene_validity` |
| `ClinVar_get_variant_details` | `dnakit.databases.clin_var_get_variant_details` |
| `ClinVar_search_variants` | `dnakit.databases.clin_var_search_variants` |
| `DisGeNET_get_gda` | `dnakit.databases.dis_ge_net_get_gda` |
| `DisGeNET_get_vda` | `dnakit.databases.dis_ge_net_get_vda` |
| `DisGeNET_search_disease` | `dnakit.databases.dis_ge_net_search_disease` |
| `DisGeNET_search_gene` | `dnakit.databases.dis_ge_net_search_gene` |
| `ENCODE_get_biosample` | `dnakit.databases.encode_get_biosample` |
| `ENCODE_get_experiment` | `dnakit.databases.encode_get_experiment` |
| `ENCODE_search_experiments` | `dnakit.databases.encode_search_experiments` |
| `ESM_describe_sae_feature` | `dnakit.databases.esm_describe_sae_feature` |
| `ESM_explain_variant_mechanism` | `dnakit.databases.esm_explain_variant_mechanism` |
| `ESM_get_region_sae_features` | `dnakit.databases.esm_get_region_sae_features` |
| `ESM_score_variant_sae_batch` | `dnakit.databases.esm_score_variant_sae_batch` |
| `ESM_score_variant_sae_disruption` | `dnakit.databases.esm_score_variant_sae_disruption` |
| `EVE_get_gene_info` | `dnakit.databases.eve_get_gene_info` |
| `EVE_get_variant_score` | `dnakit.databases.eve_get_variant_score` |
| `EnsemblVEP_annotate_hgvs` | `dnakit.databases.ensembl_vep_annotate_hgvs` |
| `EnsemblVar_get_variant_consequences` | `dnakit.databases.ensembl_var_get_variant_consequences` |
| `EuropePMC_search_articles` | `dnakit.databases.europe_pmc_search_articles` |
| `Evo2_score_variant` | `dnakit.databases.evo2_score_variant` |
| `GTEx_get_median_gene_expression` | `dnakit.databases.gtex_get_median_gene_expression` |
| `InterPro_get_protein_domains` | `dnakit.databases.inter_pro_get_protein_domains` |
| `MedRxiv_get_preprint` | `dnakit.databases.med_rxiv_get_preprint` |
| `MyVariant_get_pathogenicity_scores` | `dnakit.databases.my_variant_get_pathogenicity_scores` |
| `MyVariant_query_variants` | `dnakit.databases.my_variant_query_variants` |
| `NCBIGene_search` | `dnakit.databases.ncbigene_search` |
| `NvidiaNIM_alphafold2` | `dnakit.databases.nvidia_nim_alphafold2` |
| `OMIM_get_clinical_synopsis` | `dnakit.databases.omim_get_clinical_synopsis` |
| `OMIM_get_entry` | `dnakit.databases.omim_get_entry` |
| `OMIM_get_gene_map` | `dnakit.databases.omim_get_gene_map` |
| `OMIM_search` | `dnakit.databases.omim_search` |
| `PDBe_get_uniprot_mappings` | `dnakit.databases.pdbe_get_uniprot_mappings` |
| `PubMed_get_article` | `dnakit.databases.pub_med_get_article` |
| `PubMed_search_articles` | `dnakit.databases.pub_med_search_articles` |
| `RCSBData_get_entry` | `dnakit.databases.rcsbdata_get_entry` |
| `SemanticScholar_search_papers` | `dnakit.databases.semantic_scholar_search_papers` |
| `SpliceAI_get_max_delta` | `dnakit.databases.splice_ai_get_max_delta` |
| `SpliceAI_predict_pangolin` | `dnakit.databases.splice_ai_predict_pangolin` |
| `SpliceAI_predict_splice` | `dnakit.databases.splice_ai_predict_splice` |
| `Tark_get_mane_transcripts` | `dnakit.databases.tark_get_mane_transcripts` |
| `Tark_get_transcript` | `dnakit.databases.tark_get_transcript` |
| `UniProt_get_function_by_accession` | `dnakit.databases.uni_prot_get_function_by_accession` |
| `UniProt_get_sequence_by_accession` | `dnakit.databases.uni_prot_get_sequence_by_accession` |
| `VariantValidator_gene2transcripts` | `dnakit.databases.variant_validator_gene2transcripts` |
| `VariantValidator_validate_variant` | `dnakit.databases.variant_validator_validate_variant` |
| `alphafold_get_prediction` | `dnakit.databases.alphafold_get_prediction` |
| `civic_get_variants_by_gene` | `dnakit.databases.civic_get_variants_by_gene` |
| `civic_search_assertions` | `dnakit.databases.civic_search_assertions` |
| `civic_search_evidence_items` | `dnakit.databases.civic_search_evidence_items` |
| `civic_search_genes` | `dnakit.databases.civic_search_genes` |
| `ensembl_lookup_gene` | `dnakit.databases.ensembl_lookup_gene` |
| `gnomad_get_variant` | `dnakit.databases.gnomad_get_variant` |
| `gnomad_search_variants` | `dnakit.databases.gnomad_search_variants` |
| `openalex_search_works` | `dnakit.databases.openalex_search_works` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.

## Regulatory result types

`ch_ipatlas_get_peak_data(gene=...)` resolves the GRCh38 gene locus and obtains an assembled track link. A link response has `status="requires_data"`; it contains no observed peaks. Download and verify the selected BED, then pass `peak_file=...` to compute interval overlaps. `flank_bp` and `limit` control the queried interval and returned rows. This is not an antigen-name search.

The inherited `ch_ipatlas_enrichment_analysis` adapter supplies web-submission instructions, not computed enrichment scores; DNAKit marks this `requires_submission`. Complete that external analysis and inspect actual statistics before interpreting TF binding. Peak overlap alone does not establish motif disruption. See [ChIP-Atlas's API](https://chip-atlas.org/agents).
