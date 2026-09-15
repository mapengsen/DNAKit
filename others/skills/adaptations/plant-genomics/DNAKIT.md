# DNAKit integration

Canonical skill: `plant-genomics`.

```bash
python -m pip install dnakit
dnakit skills doctor plant-genomics
dnakit skills install plant-genomics --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('ensembl_compara_get_orthologues')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `EnsemblCompara_get_orthologues` | `dnakit.databases.ensembl_compara_get_orthologues` |
| `KEGG_get_pathway_genes` | `dnakit.databases.kegg_get_pathway_genes` |
| `POWO_search_plants` | `dnakit.databases.powo_search_plants` |
| `PlantReactome_get_pathway` | `dnakit.databases.plant_reactome_get_pathway` |
| `PlantReactome_list_species` | `dnakit.databases.plant_reactome_list_species` |
| `PlantReactome_search_pathways` | `dnakit.databases.plant_reactome_search_pathways` |
| `PubMed_search_articles` | `dnakit.databases.pub_med_search_articles` |
| `USDA_plants_get_characteristics` | `dnakit.databases.usda_plants_get_characteristics` |
| `USDA_plants_get_profile` | `dnakit.databases.usda_plants_get_profile` |
| `UniProt_get_function_by_accession` | `dnakit.databases.uni_prot_get_function_by_accession` |
| `UniProt_search` | `dnakit.databases.uni_prot_search` |
| `ensembl_lookup_gene` | `dnakit.databases.ensembl_lookup_gene` |
| `kegg_search_pathway` | `dnakit.databases.kegg_search_pathway` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
