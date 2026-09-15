# DNAKit integration

Canonical skill: `gene-regulatory-networks`.

```bash
python -m pip install dnakit
dnakit skills doctor gene-regulatory-networks
dnakit skills install gene-regulatory-networks --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('bio_grid_get_interactions')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `BioGRID_get_interactions` | `dnakit.databases.bio_grid_get_interactions` |
| `ENCODE_search_histone_experiments` | `dnakit.databases.encode_search_histone_experiments` |
| `EuropePMC_search_articles` | `dnakit.databases.europe_pmc_search_articles` |
| `GTEx_query_eqtl` | `dnakit.databases.gtex_query_eqtl` |
| `PubMed_search_articles` | `dnakit.databases.pub_med_search_articles` |
| `RegulomeDB_query_variant` | `dnakit.databases.regulome_db_query_variant` |
| `STRING_functional_enrichment` | `dnakit.databases.string_functional_enrichment` |
| `STRING_get_interaction_partners` | `dnakit.databases.string_get_interaction_partners` |
| `STRING_get_network` | `dnakit.databases.string_get_network` |
| `enrichr_gene_enrichment_analysis` | `dnakit.databases.enrichr_gene_enrichment_analysis` |
| `intact_get_interaction_network` | `dnakit.databases.intact_get_interaction_network` |
| `jaspar_get_matrix` | `dnakit.databases.jaspar_get_matrix` |
| `jaspar_search_matrices` | `dnakit.databases.jaspar_search_matrices` |
| `ols_search_terms` | `dnakit.databases.ols_search_terms` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
