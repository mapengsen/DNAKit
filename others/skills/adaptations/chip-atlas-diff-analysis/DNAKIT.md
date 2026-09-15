# DNAKit integration

Canonical skill: `chip-atlas-diff-analysis`.

```bash
python -m pip install dnakit
dnakit skills doctor chip-atlas-diff-analysis
dnakit skills install chip-atlas-diff-analysis --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts chip-atlas-diff-analysis
dnakit skills run chip-atlas-diff-analysis --script scripts/__init__.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `scripts/__init__.py`
- `scripts/annotate_genes.py`
- `scripts/export_all.py`
- `scripts/filter_regions.py`
- `scripts/generate_all_plots.py`
- `scripts/load_example_data.py`
- `scripts/load_user_data.py`
- `scripts/parse_diff_bed.py`
- `scripts/qc_checks.py`
- `scripts/query_chip_atlas_diff.py`
- `scripts/run_diff_workflow.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
