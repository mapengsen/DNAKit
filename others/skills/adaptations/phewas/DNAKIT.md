# DNAKit integration

Canonical skill: `phewas`.

```bash
python -m pip install dnakit
dnakit skills doctor phewas
dnakit skills install phewas --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('bio_bank_japan_phewas_by_variant')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `BioBankJapan_phewas_by_variant` | `dnakit.databases.bio_bank_japan_phewas_by_variant` |
| `FinnGen_get_variant_finemapping` | `dnakit.databases.finn_gen_get_variant_finemapping` |
| `Genebass_gene_burden_phewas` | `dnakit.databases.genebass_gene_burden_phewas` |
| `TPMI_phewas_by_variant` | `dnakit.databases.tpmi_phewas_by_variant` |
| `UKBTOPMed_phewas_by_variant` | `dnakit.databases.ukbtopmed_phewas_by_variant` |

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
