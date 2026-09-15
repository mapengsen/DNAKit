# DNAKit integration

Canonical skill: `gwas-trait-to-gene`.

```bash
python -m pip install dnakit
dnakit skills doctor gwas-trait-to-gene
dnakit skills install gwas-trait-to-gene --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('open_targets_get_credible_set_detail')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `OpenTargets_get_credible_set_detail` | `dnakit.databases.open_targets_get_credible_set_detail` |
| `OpenTargets_get_gwas_study` | `dnakit.databases.open_targets_get_gwas_study` |
| `OpenTargets_get_study_credible_sets` | `dnakit.databases.open_targets_get_study_credible_sets` |
| `OpenTargets_get_variant_credible_sets` | `dnakit.databases.open_targets_get_variant_credible_sets` |
| `OpenTargets_get_variant_info` | `dnakit.databases.open_targets_get_variant_info` |
| `OpenTargets_search_gwas_studies_by_disease` | `dnakit.databases.open_targets_search_gwas_studies_by_disease` |
| `gwas_get_associations_for_snp` | `dnakit.databases.gwas_get_associations_for_snp` |
| `gwas_get_associations_for_study` | `dnakit.databases.gwas_get_associations_for_study` |
| `gwas_get_associations_for_trait` | `dnakit.databases.gwas_get_associations_for_trait` |
| `gwas_get_snp_by_id` | `dnakit.databases.gwas_get_snp_by_id` |
| `gwas_get_snps_for_gene` | `dnakit.databases.gwas_get_snps_for_gene` |
| `gwas_get_studies_for_trait` | `dnakit.databases.gwas_get_studies_for_trait` |
| `gwas_get_study_by_id` | `dnakit.databases.gwas_get_study_by_id` |
| `gwas_get_variants_for_trait` | `dnakit.databases.gwas_get_variants_for_trait` |
| `gwas_search_associations` | `dnakit.databases.gwas_search_associations` |
| `gwas_search_snps` | `dnakit.databases.gwas_search_snps` |
| `gwas_search_studies` | `dnakit.databases.gwas_search_studies` |

## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts gwas-trait-to-gene
dnakit skills run gwas-trait-to-gene --script python_implementation.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `python_implementation.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
