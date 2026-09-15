# DNAKit integration

Canonical skill: `microbial-genome-characterization`.

```bash
python -m pip install dnakit
dnakit skills doctor microbial-genome-characterization
dnakit skills install microbial-genome-characterization --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('ncbidatasets_get_genome_assembly')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `NCBIDatasets_get_genome_assembly` | `dnakit.databases.ncbidatasets_get_genome_assembly` |
| `NCBIDatasets_get_sequence_reports` | `dnakit.databases.ncbidatasets_get_sequence_reports` |
| `NCBIDatasets_get_taxonomy` | `dnakit.databases.ncbidatasets_get_taxonomy` |
| `NCBIDatasets_list_genomes_by_taxon` | `dnakit.databases.ncbidatasets_list_genomes_by_taxon` |
| `NCBIDatasets_suggest_taxonomy` | `dnakit.databases.ncbidatasets_suggest_taxonomy` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
