# DNAKit integration

Canonical skill: `comparative-genomics`.

```bash
python -m pip install dnakit
dnakit skills doctor comparative-genomics
dnakit skills install comparative-genomics --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('blast_protein_search')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `BLAST_protein_search` | `dnakit.databases.blast_protein_search` |
| `EnsemblCompara_get_gene_tree` | `dnakit.databases.ensembl_compara_get_gene_tree` |
| `EnsemblCompara_get_orthologues` | `dnakit.databases.ensembl_compara_get_orthologues` |
| `EnsemblCompara_get_paralogues` | `dnakit.databases.ensembl_compara_get_paralogues` |
| `EnsemblSeq_get_id_sequence` | `dnakit.databases.ensembl_seq_get_id_sequence` |
| `MonarchV3_get_associations` | `dnakit.databases.monarch_v3_get_associations` |
| `Monarch_get_gene_diseases` | `dnakit.databases.monarch_get_gene_diseases` |
| `Monarch_get_gene_phenotypes` | `dnakit.databases.monarch_get_gene_phenotypes` |
| `Monarch_search_gene` | `dnakit.databases.monarch_search_gene` |
| `NCBIDatasets_get_genome_assembly` | `dnakit.databases.ncbidatasets_get_genome_assembly` |
| `NCBIDatasets_get_sequence_reports` | `dnakit.databases.ncbidatasets_get_sequence_reports` |
| `NCBIDatasets_list_genomes_by_taxon` | `dnakit.databases.ncbidatasets_list_genomes_by_taxon` |
| `NCBI_fetch_accessions` | `dnakit.databases.ncbi_fetch_accessions` |
| `NCBI_get_sequence` | `dnakit.databases.ncbi_get_sequence` |
| `NCBI_search_nucleotide` | `dnakit.databases.ncbi_search_nucleotide` |
| `OpenTargets_get_target_homologues_by_ensemblID` | `dnakit.databases.open_targets_get_target_homologues_by_ensembl_id` |
| `Sequence_dn_ds` | `dnakit.databases.sequence_dn_ds` |
| `UniProt_get_function_by_accession` | `dnakit.databases.uni_prot_get_function_by_accession` |
| `UniProt_search` | `dnakit.databases.uni_prot_search` |
| `ensembl_get_homology` | `dnakit.databases.ensembl_get_homology` |
| `ensembl_lookup_gene` | `dnakit.databases.ensembl_lookup_gene` |

## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts comparative-genomics
dnakit skills run comparative-genomics --script scripts/dnds.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `scripts/dnds.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
