# DNAKit integration

Canonical skill: `epigenomics`.

```bash
python -m pip install dnakit
dnakit skills doctor epigenomics
dnakit skills install epigenomics --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('ch_ipatlas_enrichment_analysis')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `ChIPAtlas_enrichment_analysis` | `dnakit.databases.ch_ipatlas_enrichment_analysis` |
| `ChIPAtlas_get_experiments` | `dnakit.databases.ch_ipatlas_get_experiments` |
| `ChIPAtlas_get_peak_data` | `dnakit.databases.ch_ipatlas_get_peak_data` |
| `ChIPAtlas_search_datasets` | `dnakit.databases.ch_ipatlas_search_datasets` |
| `ENCODE_search_experiments` | `dnakit.databases.encode_search_experiments` |
| `ENCODE_search_histone_experiments` | `dnakit.databases.encode_search_histone_experiments` |
| `ENCODE_search_rnaseq_experiments` | `dnakit.databases.encode_search_rnaseq_experiments` |
| `FourDN_search_data` | `dnakit.databases.four_dn_search_data` |
| `GEO_search_atacseq_datasets` | `dnakit.databases.geo_search_atacseq_datasets` |
| `GEO_search_rnaseq_datasets` | `dnakit.databases.geo_search_rnaseq_datasets` |
| `GO_get_annotations_for_gene` | `dnakit.databases.go_get_annotations_for_gene` |
| `GTEx_get_median_gene_expression` | `dnakit.databases.gtex_get_median_gene_expression` |
| `GTEx_query_eqtl` | `dnakit.databases.gtex_query_eqtl` |
| `HGNC_fetch_gene_by_symbol` | `dnakit.databases.hgnc_fetch_gene_by_symbol` |
| `MyGene_batch_query` | `dnakit.databases.my_gene_batch_query` |
| `MyGene_query_genes` | `dnakit.databases.my_gene_query_genes` |
| `ReMap_get_transcription_factor_binding` | `dnakit.databases.re_map_get_transcription_factor_binding` |
| `RegulomeDB_query_variant` | `dnakit.databases.regulome_db_query_variant` |
| `SCREEN_get_regulatory_elements` | `dnakit.databases.screen_get_regulatory_elements` |
| `SRA_get_experiment` | `dnakit.databases.sra_get_experiment` |
| `SRA_search_experiments` | `dnakit.databases.sra_search_experiments` |
| `ensembl_get_overlap_features` | `dnakit.databases.ensembl_get_overlap_features` |
| `ensembl_get_regulatory_features` | `dnakit.databases.ensembl_get_regulatory_features` |
| `ensembl_lookup_gene` | `dnakit.databases.ensembl_lookup_gene` |
| `jaspar_search_matrices` | `dnakit.databases.jaspar_search_matrices` |
| `read_executed_notebook` | `dnakit.databases.read_executed_notebook` |

## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts epigenomics
dnakit skills run epigenomics --script scripts/methylation_density.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `scripts/methylation_density.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.

## Regulatory result types

`ch_ipatlas_get_peak_data(gene=...)` resolves the GRCh38 gene locus and obtains an assembled track link. A link response has `status="requires_data"`; it contains no observed peaks. Download and verify the selected BED, then pass `peak_file=...` to compute interval overlaps. `flank_bp` and `limit` control the queried interval and returned rows. This is not an antigen-name search.

The inherited `ch_ipatlas_enrichment_analysis` adapter supplies web-submission instructions, not computed enrichment scores; DNAKit marks this `requires_submission`. Complete that external analysis and inspect actual statistics before interpreting TF binding. Peak overlap alone does not establish motif disruption. See [ChIP-Atlas's API](https://chip-atlas.org/agents).
