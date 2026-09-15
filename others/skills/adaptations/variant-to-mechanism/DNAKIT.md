# DNAKit integration

Canonical skill: `variant-to-mechanism`.

```bash
python -m pip install dnakit
dnakit skills doctor variant-to-mechanism
dnakit skills install variant-to-mechanism --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('dis_ge_net_search_gene')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `DisGeNET_search_gene` | `dnakit.databases.dis_ge_net_search_gene` |
| `ENCODE_search_chromatin_accessibility` | `dnakit.databases.encode_search_chromatin_accessibility` |
| `ENCODE_search_histone_experiments` | `dnakit.databases.encode_search_histone_experiments` |
| `EnsemblVEP_annotate_rsid` | `dnakit.databases.ensembl_vep_annotate_rsid` |
| `FAVOR_annotate_variant` | `dnakit.databases.favor_annotate_variant` |
| `GTEx_get_median_gene_expression` | `dnakit.databases.gtex_get_median_gene_expression` |
| `GTEx_query_eqtl` | `dnakit.databases.gtex_query_eqtl` |
| `GenCC_search_gene` | `dnakit.databases.gen_cc_search_gene` |
| `MyGene_query_genes` | `dnakit.databases.my_gene_query_genes` |
| `MyVariant_query_variants` | `dnakit.databases.my_variant_query_variants` |
| `OpenTargets_get_diseases_phenotypes_by_target_ensembl` | `dnakit.databases.open_targets_get_diseases_phenotypes_by_target_ensembl` |
| `OpenTargets_get_variant_credible_sets` | `dnakit.databases.open_targets_get_variant_credible_sets` |
| `OpenTargets_search_gwas_studies_by_disease` | `dnakit.databases.open_targets_search_gwas_studies_by_disease` |
| `OpenTargets_target_disease_evidence` | `dnakit.databases.open_targets_target_disease_evidence` |
| `PANTHER_enrichment` | `dnakit.databases.panther_enrichment` |
| `ReactomeAnalysis_pathway_enrichment` | `dnakit.databases.reactome_analysis_pathway_enrichment` |
| `RegulomeDB_query_variant` | `dnakit.databases.regulome_db_query_variant` |
| `STRING_get_interaction_partners` | `dnakit.databases.string_get_interaction_partners` |
| `UCSC_get_encode_cCREs` | `dnakit.databases.ucsc_get_encode_c_cres` |
| `UniProt_get_function_by_accession` | `dnakit.databases.uni_prot_get_function_by_accession` |
| `gwas_get_associations_for_trait` | `dnakit.databases.gwas_get_associations_for_trait` |
| `gwas_search_associations` | `dnakit.databases.gwas_search_associations` |
| `gwas_search_snps` | `dnakit.databases.gwas_search_snps` |
| `ols_search_terms` | `dnakit.databases.ols_search_terms` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
