# DNAKit integration

Canonical skill: `regulatory-variant-analysis`.

```bash
python -m pip install dnakit
dnakit skills doctor regulatory-variant-analysis
dnakit skills install regulatory-variant-analysis --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('encode_search_histone_experiments')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `ENCODE_search_histone_experiments` | `dnakit.databases.encode_search_histone_experiments` |
| `EnsemblVEP_annotate_rsid` | `dnakit.databases.ensembl_vep_annotate_rsid` |
| `FAVOR_annotate_variant` | `dnakit.databases.favor_annotate_variant` |
| `GTEx_get_median_gene_expression` | `dnakit.databases.gtex_get_median_gene_expression` |
| `GTEx_query_eqtl` | `dnakit.databases.gtex_query_eqtl` |
| `OpenTargets_get_disease_id_description_by_name` | `dnakit.databases.open_targets_get_disease_id_description_by_name` |
| `OpenTargets_multi_entity_search_by_query_string` | `dnakit.databases.open_targets_multi_entity_search_by_query_string` |
| `OpenTargets_search_gwas_studies_by_disease` | `dnakit.databases.open_targets_search_gwas_studies_by_disease` |
| `RegulomeDB_query_variant` | `dnakit.databases.regulome_db_query_variant` |
| `gwas_get_snps_for_gene` | `dnakit.databases.gwas_get_snps_for_gene` |
| `gwas_get_variants_for_trait` | `dnakit.databases.gwas_get_variants_for_trait` |
| `gwas_search_associations` | `dnakit.databases.gwas_search_associations` |
| `ols_search_terms` | `dnakit.databases.ols_search_terms` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
