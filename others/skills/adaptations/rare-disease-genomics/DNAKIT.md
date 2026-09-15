# DNAKit integration

Canonical skill: `rare-disease-genomics`.

```bash
python -m pip install dnakit
dnakit skills doctor rare-disease-genomics
dnakit skills install rare-disease-genomics --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('clin_var_get_clinical_significance')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `ClinVar_get_clinical_significance` | `dnakit.databases.clin_var_get_clinical_significance` |
| `ClinVar_get_variant_details` | `dnakit.databases.clin_var_get_variant_details` |
| `ClinVar_search_variants` | `dnakit.databases.clin_var_search_variants` |
| `EuropePMC_search_articles` | `dnakit.databases.europe_pmc_search_articles` |
| `GenCC_search_disease` | `dnakit.databases.gen_cc_search_disease` |
| `GenCC_search_gene` | `dnakit.databases.gen_cc_search_gene` |
| `HMDB_get_diseases` | `dnakit.databases.hmdb_get_diseases` |
| `HMDB_get_metabolite` | `dnakit.databases.hmdb_get_metabolite` |
| `HMDB_search` | `dnakit.databases.hmdb_search` |
| `Orphanet_get_classification` | `dnakit.databases.orphanet_get_classification` |
| `Orphanet_get_disease` | `dnakit.databases.orphanet_get_disease` |
| `Orphanet_get_epidemiology` | `dnakit.databases.orphanet_get_epidemiology` |
| `Orphanet_get_gene_diseases` | `dnakit.databases.orphanet_get_gene_diseases` |
| `Orphanet_get_genes` | `dnakit.databases.orphanet_get_genes` |
| `Orphanet_get_icd_mapping` | `dnakit.databases.orphanet_get_icd_mapping` |
| `Orphanet_get_natural_history` | `dnakit.databases.orphanet_get_natural_history` |
| `Orphanet_get_phenotypes` | `dnakit.databases.orphanet_get_phenotypes` |
| `Orphanet_search_diseases` | `dnakit.databases.orphanet_search_diseases` |
| `ols_get_term_ancestors` | `dnakit.databases.ols_get_term_ancestors` |
| `ols_get_term_children` | `dnakit.databases.ols_get_term_children` |
| `ols_get_term_info` | `dnakit.databases.ols_get_term_info` |
| `ols_search_terms` | `dnakit.databases.ols_search_terms` |
| `search_clinical_trials` | `dnakit.databases.search_clinical_trials` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
