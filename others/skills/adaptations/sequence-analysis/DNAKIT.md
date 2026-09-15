# DNAKit integration

Canonical skill: `sequence-analysis`.

```bash
python -m pip install dnakit
dnakit skills doctor sequence-analysis
dnakit skills install sequence-analysis --agent codex
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
| `ClinVar_search_variants` | `dnakit.databases.clin_var_search_variants` |
| `DNA_translate_reading_frames` | `dnakit.databases.dna_translate_reading_frames` |
| `EnsemblCompara_get_orthologues` | `dnakit.databases.ensembl_compara_get_orthologues` |
| `EnsemblSeq_get_id_sequence` | `dnakit.databases.ensembl_seq_get_id_sequence` |
| `EnsemblSeq_get_region_sequence` | `dnakit.databases.ensembl_seq_get_region_sequence` |
| `EnsemblVEP_annotate_hgvs` | `dnakit.databases.ensembl_vep_annotate_hgvs` |
| `InterPro_get_entries_for_protein` | `dnakit.databases.inter_pro_get_entries_for_protein` |
| `NCBIDatasets_get_gene` | `dnakit.databases.ncbidatasets_get_gene` |
| `NCBIDatasets_get_gene_by_symbol` | `dnakit.databases.ncbidatasets_get_gene_by_symbol` |
| `NCBIDatasets_get_orthologs` | `dnakit.databases.ncbidatasets_get_orthologs` |
| `NCBIGene_get_summary` | `dnakit.databases.ncbigene_get_summary` |
| `NCBIGene_search` | `dnakit.databases.ncbigene_search` |
| `NCBIProtein_get_summary` | `dnakit.databases.ncbiprotein_get_summary` |
| `NCBI_fetch_accessions` | `dnakit.databases.ncbi_fetch_accessions` |
| `NCBI_get_sequence` | `dnakit.databases.ncbi_get_sequence` |
| `NCBI_search_nucleotide` | `dnakit.databases.ncbi_search_nucleotide` |
| `Pfam_get_protein_annotations` | `dnakit.databases.pfam_get_protein_annotations` |
| `PubMed_search_articles` | `dnakit.databases.pub_med_search_articles` |
| `Sequence_count_residues` | `dnakit.databases.sequence_count_residues` |
| `Sequence_gc_content` | `dnakit.databases.sequence_gc_content` |
| `Sequence_reverse_complement` | `dnakit.databases.sequence_reverse_complement` |
| `Sequence_stats` | `dnakit.databases.sequence_stats` |
| `UniProt_get_entry_by_accession` | `dnakit.databases.uni_prot_get_entry_by_accession` |
| `UniProt_get_sequence_by_accession` | `dnakit.databases.uni_prot_get_sequence_by_accession` |
| `UniProt_search` | `dnakit.databases.uni_prot_search` |
| `ensembl_get_sequence` | `dnakit.databases.ensembl_get_sequence` |
| `proteins_api_get_features` | `dnakit.databases.proteins_api_get_features` |
| `proteins_api_search` | `dnakit.databases.proteins_api_search` |
| `read_executed_notebook` | `dnakit.databases.read_executed_notebook` |

## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts sequence-analysis
dnakit skills run sequence-analysis --script scripts/amino_acids.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `scripts/amino_acids.py`
- `scripts/biology_facts.py`
- `scripts/sequence_tools.py`
- `scripts/translate_dna.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
