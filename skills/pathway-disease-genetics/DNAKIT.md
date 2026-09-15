# DNAKit integration

Canonical skill: `pathway-disease-genetics`.

```bash
python -m pip install dnakit
dnakit skills doctor pathway-disease-genetics
dnakit skills install pathway-disease-genetics --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('dgidb_get_drug_gene_interactions')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `DGIdb_get_drug_gene_interactions` | `dnakit.databases.dgidb_get_drug_gene_interactions` |
| `DGIdb_get_gene_druggability` | `dnakit.databases.dgidb_get_gene_druggability` |
| `EnsemblVEP_annotate_rsid` | `dnakit.databases.ensembl_vep_annotate_rsid` |
| `GTEx_get_expression_summary` | `dnakit.databases.gtex_get_expression_summary` |
| `GTEx_get_median_gene_expression` | `dnakit.databases.gtex_get_median_gene_expression` |
| `GTEx_query_eqtl` | `dnakit.databases.gtex_query_eqtl` |
| `KEGG_get_gene_pathways` | `dnakit.databases.kegg_get_gene_pathways` |
| `KEGG_get_pathway_genes` | `dnakit.databases.kegg_get_pathway_genes` |
| `OpenTargets_get_associated_drugs_by_target_ensemblID` | `dnakit.databases.open_targets_get_associated_drugs_by_target_ensembl_id` |
| `OpenTargets_multi_entity_search_by_query_string` | `dnakit.databases.open_targets_multi_entity_search_by_query_string` |
| `OpenTargets_target_disease_evidence` | `dnakit.databases.open_targets_target_disease_evidence` |
| `PANTHER_enrichment` | `dnakit.databases.panther_enrichment` |
| `ReactomeAnalysis_pathway_enrichment` | `dnakit.databases.reactome_analysis_pathway_enrichment` |
| `Reactome_get_participants` | `dnakit.databases.reactome_get_participants` |
| `Reactome_map_uniprot_to_pathways` | `dnakit.databases.reactome_map_uniprot_to_pathways` |
| `STRING_functional_enrichment` | `dnakit.databases.string_functional_enrichment` |
| `gwas_get_associations_for_trait` | `dnakit.databases.gwas_get_associations_for_trait` |
| `gwas_get_snps_for_gene` | `dnakit.databases.gwas_get_snps_for_gene` |
| `gwas_get_variants_for_trait` | `dnakit.databases.gwas_get_variants_for_trait` |
| `gwas_search_associations` | `dnakit.databases.gwas_search_associations` |
| `humanbase_ppi_analysis` | `dnakit.databases.humanbase_ppi_analysis` |
| `kegg_find_genes` | `dnakit.databases.kegg_find_genes` |
| `kegg_search_pathway` | `dnakit.databases.kegg_search_pathway` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
