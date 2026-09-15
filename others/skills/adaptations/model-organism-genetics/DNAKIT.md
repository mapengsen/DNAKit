# DNAKit integration

Canonical skill: `model-organism-genetics`.

```bash
python -m pip install dnakit
dnakit skills doctor model-organism-genetics
dnakit skills install model-organism-genetics --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('clin_gen_search_gene_validity')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `ClinGen_search_gene_validity` | `dnakit.databases.clin_gen_search_gene_validity` |
| `ClinVar_search_variants` | `dnakit.databases.clin_var_search_variants` |
| `EnsemblCompara_get_orthologues` | `dnakit.databases.ensembl_compara_get_orthologues` |
| `FlyBase_get_gene` | `dnakit.databases.fly_base_get_gene` |
| `FlyBase_get_gene_alleles` | `dnakit.databases.fly_base_get_gene_alleles` |
| `FlyBase_get_gene_disease_models` | `dnakit.databases.fly_base_get_gene_disease_models` |
| `FlyBase_get_gene_expression` | `dnakit.databases.fly_base_get_gene_expression` |
| `FlyBase_get_gene_interactions` | `dnakit.databases.fly_base_get_gene_interactions` |
| `FlyBase_get_gene_orthologs` | `dnakit.databases.fly_base_get_gene_orthologs` |
| `FlyMine_search` | `dnakit.databases.fly_mine_search` |
| `GBIF_get_species` | `dnakit.databases.gbif_get_species` |
| `GBIF_search_species` | `dnakit.databases.gbif_search_species` |
| `HPO_search_terms` | `dnakit.databases.hpo_search_terms` |
| `MGI_get_gene` | `dnakit.databases.mgi_get_gene` |
| `MGI_get_phenotypes` | `dnakit.databases.mgi_get_phenotypes` |
| `MGI_search_genes` | `dnakit.databases.mgi_search_genes` |
| `MonarchV3_get_associations` | `dnakit.databases.monarch_v3_get_associations` |
| `MonarchV3_phenotype_similarity_search` | `dnakit.databases.monarch_v3_phenotype_similarity_search` |
| `Monarch_search_gene` | `dnakit.databases.monarch_search_gene` |
| `MyGene_query_genes` | `dnakit.databases.my_gene_query_genes` |
| `NCBIDatasets_get_orthologs` | `dnakit.databases.ncbidatasets_get_orthologs` |
| `OMIM_search` | `dnakit.databases.omim_search` |
| `PANTHER_ortholog` | `dnakit.databases.panther_ortholog` |
| `RGD_get_annotations` | `dnakit.databases.rgd_get_annotations` |
| `RGD_get_gene` | `dnakit.databases.rgd_get_gene` |
| `RGD_get_orthologs` | `dnakit.databases.rgd_get_orthologs` |
| `RGD_get_qtls_in_region` | `dnakit.databases.rgd_get_qtls_in_region` |
| `RGD_get_strain` | `dnakit.databases.rgd_get_strain` |
| `RGD_get_strain_annotations` | `dnakit.databases.rgd_get_strain_annotations` |
| `RGD_search_genes` | `dnakit.databases.rgd_search_genes` |
| `RGD_search_strains` | `dnakit.databases.rgd_search_strains` |
| `ReactomeAnalysis_pathway_enrichment` | `dnakit.databases.reactome_analysis_pathway_enrichment` |
| `SGD_get_gene` | `dnakit.databases.sgd_get_gene` |
| `SGD_get_go_annotations` | `dnakit.databases.sgd_get_go_annotations` |
| `SGD_get_interactions` | `dnakit.databases.sgd_get_interactions` |
| `SGD_get_phenotypes` | `dnakit.databases.sgd_get_phenotypes` |
| `SGD_search` | `dnakit.databases.sgd_search` |
| `STRING_get_network` | `dnakit.databases.string_get_network` |
| `UniProt_search` | `dnakit.databases.uni_prot_search` |
| `WormBase_get_expression` | `dnakit.databases.worm_base_get_expression` |
| `WormBase_get_gene` | `dnakit.databases.worm_base_get_gene` |
| `WormBase_get_phenotypes` | `dnakit.databases.worm_base_get_phenotypes` |
| `Xenbase_get_gene` | `dnakit.databases.xenbase_get_gene` |
| `Xenbase_search_genes` | `dnakit.databases.xenbase_search_genes` |
| `ZFIN_get_gene` | `dnakit.databases.zfin_get_gene` |
| `ZFIN_get_gene_expression` | `dnakit.databases.zfin_get_gene_expression` |
| `ZFIN_get_gene_phenotypes` | `dnakit.databases.zfin_get_gene_phenotypes` |
| `ensembl_lookup_gene` | `dnakit.databases.ensembl_lookup_gene` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
