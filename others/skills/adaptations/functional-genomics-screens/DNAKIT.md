# DNAKit integration

Canonical skill: `functional-genomics-screens`.

```bash
python -m pip install dnakit
dnakit skills doctor functional-genomics-screens
dnakit skills install functional-genomics-screens --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('cosmic_get_mutations_by_gene')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `COSMIC_get_mutations_by_gene` | `dnakit.databases.cosmic_get_mutations_by_gene` |
| `DGIdb_get_drug_gene_interactions` | `dnakit.databases.dgidb_get_drug_gene_interactions` |
| `DGIdb_get_gene_druggability` | `dnakit.databases.dgidb_get_gene_druggability` |
| `DepMap_get_gene_dependencies` | `dnakit.databases.dep_map_get_gene_dependencies` |
| `DepMap_search_cell_lines` | `dnakit.databases.dep_map_search_cell_lines` |
| `PubMed_search_articles` | `dnakit.databases.pub_med_search_articles` |
| `ReactomeAnalysis_pathway_enrichment` | `dnakit.databases.reactome_analysis_pathway_enrichment` |
| `STRING_functional_enrichment` | `dnakit.databases.string_functional_enrichment` |
| `STRING_get_network` | `dnakit.databases.string_get_network` |
| `UniProt_get_function_by_accession` | `dnakit.databases.uni_prot_get_function_by_accession` |
| `civic_search_evidence_items` | `dnakit.databases.civic_search_evidence_items` |
| `gnomad_get_gene_constraints` | `dnakit.databases.gnomad_get_gene_constraints` |
| `search_clinical_trials` | `dnakit.databases.search_clinical_trials` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
