# DNAKit integration

Canonical skill: `regulatory-genomics`.

```bash
python -m pip install dnakit
dnakit skills doctor regulatory-genomics
dnakit skills install regulatory-genomics --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('alpha_genome_predict_interval')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `AlphaGenome_predict_interval` | `dnakit.databases.alpha_genome_predict_interval` |
| `AlphaGenome_score_variant` | `dnakit.databases.alpha_genome_score_variant` |
| `ENCODE_get_chromatin_state` | `dnakit.databases.encode_get_chromatin_state` |
| `ENCODE_get_experiment` | `dnakit.databases.encode_get_experiment` |
| `ENCODE_search_annotations` | `dnakit.databases.encode_search_annotations` |
| `ENCODE_search_biosamples` | `dnakit.databases.encode_search_biosamples` |
| `ENCODE_search_chromatin_accessibility` | `dnakit.databases.encode_search_chromatin_accessibility` |
| `ENCODE_search_experiments` | `dnakit.databases.encode_search_experiments` |
| `ENCODE_search_histone_experiments` | `dnakit.databases.encode_search_histone_experiments` |
| `Evo2_score_variant` | `dnakit.databases.evo2_score_variant` |
| `JASPAR_get_transcription_factors` | `dnakit.databases.jaspar_get_transcription_factors` |
| `MyGene_query_genes` | `dnakit.databases.my_gene_query_genes` |
| `RegulomeDB_query_variant` | `dnakit.databases.regulome_db_query_variant` |
| `UCSC_get_encode_cCREs` | `dnakit.databases.ucsc_get_encode_c_cres` |
| `UniBind_get_dataset` | `dnakit.databases.uni_bind_get_dataset` |
| `UniBind_list_tfs` | `dnakit.databases.uni_bind_list_tfs` |
| `UniBind_search_datasets` | `dnakit.databases.uni_bind_search_datasets` |
| `ensembl_lookup_gene` | `dnakit.databases.ensembl_lookup_gene` |
| `jaspar_get_matrix` | `dnakit.databases.jaspar_get_matrix` |
| `jaspar_get_matrix_versions` | `dnakit.databases.jaspar_get_matrix_versions` |
| `jaspar_search_matrices` | `dnakit.databases.jaspar_search_matrices` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
