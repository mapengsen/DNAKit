# DNAKit integration

Canonical skill: `variant-analysis`.

```bash
python -m pip install dnakit
dnakit skills doctor variant-analysis
dnakit skills install variant-analysis --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('clin_gen_dosage_by_gene')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `ClinGen_dosage_by_gene` | `dnakit.databases.clin_gen_dosage_by_gene` |
| `ClinGen_dosage_region_search` | `dnakit.databases.clin_gen_dosage_region_search` |
| `EnsemblVEP_annotate_rsid` | `dnakit.databases.ensembl_vep_annotate_rsid` |
| `MyVariant_get_variant_annotation` | `dnakit.databases.my_variant_get_variant_annotation` |
| `MyVariant_query_variants` | `dnakit.databases.my_variant_query_variants` |
| `VCF_count_variants` | `dnakit.databases.vcf_count_variants` |
| `VCF_normalize` | `dnakit.databases.vcf_normalize` |
| `VCF_summary_stats` | `dnakit.databases.vcf_summary_stats` |
| `coding_variant_fraction` | `dnakit.databases.coding_variant_fraction` |
| `dbsnp_get_variant_by_rsid` | `dnakit.databases.dbsnp_get_variant_by_rsid` |
| `ensembl_get_structural_variants` | `dnakit.databases.ensembl_get_structural_variants` |
| `ensembl_get_sv_detail` | `dnakit.databases.ensembl_get_sv_detail` |
| `ensembl_vep_region` | `dnakit.databases.ensembl_vep_region` |
| `execute_tool` | `dnakit.databases.execute_tool` |
| `gnomad_get_sv_by_gene` | `dnakit.databases.gnomad_get_sv_by_gene` |
| `gnomad_get_sv_by_region` | `dnakit.databases.gnomad_get_sv_by_region` |
| `gnomad_get_sv_detail` | `dnakit.databases.gnomad_get_sv_detail` |
| `gnomad_get_variant` | `dnakit.databases.gnomad_get_variant` |
| `read_executed_notebook` | `dnakit.databases.read_executed_notebook` |

## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts variant-analysis
dnakit skills run variant-analysis --script python_implementation.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `python_implementation.py`
- `scripts/annotate_variants.py`
- `scripts/coding_variant_filter.py`
- `scripts/filter_variants.py`
- `scripts/gatk_haplotypecaller_pipeline.py`
- `scripts/parse_vcf.py`
- `scripts/variant_fraction.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
