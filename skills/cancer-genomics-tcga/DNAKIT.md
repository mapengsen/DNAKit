# DNAKit integration

Canonical skill: `cancer-genomics-tcga`.

```bash
python -m pip install dnakit
dnakit skills doctor cancer-genomics-tcga
dnakit skills install cancer-genomics-tcga --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('gdc_get_clinical_data')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `GDC_get_clinical_data` | `dnakit.databases.gdc_get_clinical_data` |
| `GDC_get_mutation_frequency` | `dnakit.databases.gdc_get_mutation_frequency` |
| `GDC_get_ssm_by_gene` | `dnakit.databases.gdc_get_ssm_by_gene` |
| `GDC_get_survival` | `dnakit.databases.gdc_get_survival` |
| `GDC_list_files` | `dnakit.databases.gdc_list_files` |
| `GDC_list_projects` | `dnakit.databases.gdc_list_projects` |
| `GDC_search_cases` | `dnakit.databases.gdc_search_cases` |
| `OncoKB_annotate_variant` | `dnakit.databases.onco_kb_annotate_variant` |
| `Progenetix_cnv_search` | `dnakit.databases.progenetix_cnv_search` |
| `Progenetix_list_cohorts` | `dnakit.databases.progenetix_list_cohorts` |
| `Progenetix_list_filtering_terms` | `dnakit.databases.progenetix_list_filtering_terms` |
| `Progenetix_search_biosamples` | `dnakit.databases.progenetix_search_biosamples` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
