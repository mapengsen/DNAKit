# DNAKit integration

Canonical skill: `immune-repertoire-analysis`.

```bash
python -m pip install dnakit
dnakit skills doctor immune-repertoire-analysis
dnakit skills install immune-repertoire-analysis --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('bvbrc_search_epitopes')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `BVBRC_search_epitopes` | `dnakit.databases.bvbrc_search_epitopes` |
| `ImmPort_search_studies` | `dnakit.databases.imm_port_search_studies` |
| `NCBI_SRA_search_runs` | `dnakit.databases.ncbi_sra_search_runs` |
| `PubMed_search_articles` | `dnakit.databases.pub_med_search_articles` |
| `UniProt_get_entry_by_accession` | `dnakit.databases.uni_prot_get_entry_by_accession` |
| `iedb_get_epitope_antigens` | `dnakit.databases.iedb_get_epitope_antigens` |
| `iedb_get_epitope_mhc` | `dnakit.databases.iedb_get_epitope_mhc` |
| `iedb_search_bcell` | `dnakit.databases.iedb_search_bcell` |
| `iedb_search_epitopes` | `dnakit.databases.iedb_search_epitopes` |
| `iedb_search_tcell_assays` | `dnakit.databases.iedb_search_tcell_assays` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
