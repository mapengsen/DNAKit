# DNAKit integration

Canonical skill: `hla-immunogenomics`.

```bash
python -m pip install dnakit
dnakit skills doctor hla-immunogenomics
dnakit skills install hla-immunogenomics --agent codex
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
| `DGIdb_get_drug_gene_interactions` | `dnakit.databases.dgidb_get_drug_gene_interactions` |
| `IMGT_get_gene_info` | `dnakit.databases.imgt_get_gene_info` |
| `IMGT_search_genes` | `dnakit.databases.imgt_search_genes` |
| `PubMed_search_articles` | `dnakit.databases.pub_med_search_articles` |
| `UniProt_search` | `dnakit.databases.uni_prot_search` |
| `iedb_get_epitope_mhc` | `dnakit.databases.iedb_get_epitope_mhc` |
| `iedb_search_epitopes` | `dnakit.databases.iedb_search_epitopes` |
| `iedb_search_mhc` | `dnakit.databases.iedb_search_mhc` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
