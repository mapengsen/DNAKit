# DNAKit integration

Canonical skill: `chip-atlas-peak-enrichment`.

```bash
python -m pip install dnakit
dnakit skills doctor chip-atlas-peak-enrichment
dnakit skills install chip-atlas-peak-enrichment --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts chip-atlas-peak-enrichment
dnakit skills run chip-atlas-peak-enrichment --script scripts/__init__.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `scripts/__init__.py`
- `scripts/convert_genes_to_regions.py`
- `scripts/export_all.py`
- `scripts/filter_experiments.py`
- `scripts/generate_all_plots.py`
- `scripts/load_example_data.py`
- `scripts/load_user_data.py`
- `scripts/query_chip_atlas_enrichment.py`
- `scripts/run_enrichment_workflow.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
