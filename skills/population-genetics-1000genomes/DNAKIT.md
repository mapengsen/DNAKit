# DNAKit integration

Canonical skill: `population-genetics-1000genomes`.

```bash
python -m pip install dnakit
dnakit skills doctor population-genetics-1000genomes
dnakit skills install population-genetics-1000genomes --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('igsr_list_data_collections')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `IGSR_list_data_collections` | `dnakit.databases.igsr_list_data_collections` |
| `IGSR_search_populations` | `dnakit.databases.igsr_search_populations` |
| `IGSR_search_samples` | `dnakit.databases.igsr_search_samples` |
| `gwas_get_snps_for_gene` | `dnakit.databases.gwas_get_snps_for_gene` |
| `gwas_get_variants_for_trait` | `dnakit.databases.gwas_get_variants_for_trait` |
| `gwas_search_associations` | `dnakit.databases.gwas_search_associations` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
