# DNAKit integration

Canonical skill: `microbiome-research`.

```bash
python -m pip install dnakit
dnakit skills doctor microbiome-research
dnakit skills install microbiome-research --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('ctd_get_chemical_gene_interactions')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `CTD_get_chemical_gene_interactions` | `dnakit.databases.ctd_get_chemical_gene_interactions` |
| `ENAPortal_search_samples` | `dnakit.databases.enaportal_search_samples` |
| `ENAPortal_search_studies` | `dnakit.databases.enaportal_search_studies` |
| `EuropePMC_search_articles` | `dnakit.databases.europe_pmc_search_articles` |
| `GTDB_get_species` | `dnakit.databases.gtdb_get_species` |
| `GTDB_get_taxon_info` | `dnakit.databases.gtdb_get_taxon_info` |
| `GTDB_search_genomes` | `dnakit.databases.gtdb_search_genomes` |
| `GTDB_search_taxon` | `dnakit.databases.gtdb_search_taxon` |
| `MGnify_get_genome` | `dnakit.databases.mgnify_get_genome` |
| `MGnify_get_go_terms` | `dnakit.databases.mgnify_get_go_terms` |
| `MGnify_get_interpro` | `dnakit.databases.mgnify_get_interpro` |
| `MGnify_get_study_detail` | `dnakit.databases.mgnify_get_study_detail` |
| `MGnify_get_taxonomy` | `dnakit.databases.mgnify_get_taxonomy` |
| `MGnify_list_analyses` | `dnakit.databases.mgnify_list_analyses` |
| `MGnify_list_biomes` | `dnakit.databases.mgnify_list_biomes` |
| `MGnify_search_genomes` | `dnakit.databases.mgnify_search_genomes` |
| `MGnify_search_studies` | `dnakit.databases.mgnify_search_studies` |
| `PubChem_get_CID_by_compound_name` | `dnakit.databases.pub_chem_get_cid_by_compound_name` |
| `PubChem_get_compound_properties_by_CID` | `dnakit.databases.pub_chem_get_compound_properties_by_cid` |
| `PubMed_search_articles` | `dnakit.databases.pub_med_search_articles` |
| `ReactomeAnalysis_pathway_enrichment` | `dnakit.databases.reactome_analysis_pathway_enrichment` |
| `drugbank_vocab_search` | `dnakit.databases.drugbank_vocab_search` |
| `kegg_get_pathway_info` | `dnakit.databases.kegg_get_pathway_info` |
| `kegg_search_pathway` | `dnakit.databases.kegg_search_pathway` |
| `ols_search_terms` | `dnakit.databases.ols_search_terms` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
