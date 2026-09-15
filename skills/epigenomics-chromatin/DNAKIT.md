# DNAKit integration

Canonical skill: `epigenomics-chromatin`.

```bash
python -m pip install dnakit
dnakit skills doctor epigenomics-chromatin
dnakit skills install epigenomics-chromatin --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('encode_get_chromatin_state')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `ENCODE_get_chromatin_state` | `dnakit.databases.encode_get_chromatin_state` |
| `ENCODE_search_annotations` | `dnakit.databases.encode_search_annotations` |
| `ENCODE_search_chromatin_accessibility` | `dnakit.databases.encode_search_chromatin_accessibility` |
| `ENCODE_search_histone_experiments` | `dnakit.databases.encode_search_histone_experiments` |
| `ENCODE_search_rnaseq_experiments` | `dnakit.databases.encode_search_rnaseq_experiments` |
| `GEO_search_atacseq_datasets` | `dnakit.databases.geo_search_atacseq_datasets` |
| `GEO_search_chipseq_datasets` | `dnakit.databases.geo_search_chipseq_datasets` |
| `GEO_search_rnaseq_datasets` | `dnakit.databases.geo_search_rnaseq_datasets` |
| `GTEx_calculate_eqtl` | `dnakit.databases.gtex_calculate_eqtl` |
| `GTEx_get_expression_summary` | `dnakit.databases.gtex_get_expression_summary` |
| `GTEx_get_median_gene_expression` | `dnakit.databases.gtex_get_median_gene_expression` |
| `GTEx_get_multi_tissue_eqtls` | `dnakit.databases.gtex_get_multi_tissue_eqtls` |
| `GTEx_get_single_tissue_eqtls` | `dnakit.databases.gtex_get_single_tissue_eqtls` |
| `GTEx_get_tissue_sites` | `dnakit.databases.gtex_get_tissue_sites` |
| `GTEx_query_eqtl` | `dnakit.databases.gtex_query_eqtl` |
| `MyGene_query_genes` | `dnakit.databases.my_gene_query_genes` |
| `ReMap_get_transcription_factor_binding` | `dnakit.databases.re_map_get_transcription_factor_binding` |
| `RegulomeDB_query_variant` | `dnakit.databases.regulome_db_query_variant` |
| `SCREEN_get_regulatory_elements` | `dnakit.databases.screen_get_regulatory_elements` |
| `STRING_get_functional_annotations` | `dnakit.databases.string_get_functional_annotations` |
| `UCSC_get_encode_cCREs` | `dnakit.databases.ucsc_get_encode_c_cres` |
| `eQTL_get_associations` | `dnakit.databases.e_qtl_get_associations` |
| `eQTL_list_datasets` | `dnakit.databases.e_qtl_list_datasets` |
| `ensembl_lookup_gene` | `dnakit.databases.ensembl_lookup_gene` |
| `jaspar_get_matrix` | `dnakit.databases.jaspar_get_matrix` |
| `jaspar_search_matrices` | `dnakit.databases.jaspar_search_matrices` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
