# DNAKit integration

Canonical skill: `pathogen-variant-surveillance`.

```bash
python -m pip install dnakit
dnakit skills doctor pathogen-variant-surveillance
dnakit skills install pathogen-variant-surveillance --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts pathogen-variant-surveillance
dnakit skills run pathogen-variant-surveillance --script scripts/lapis_client.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `scripts/lapis_client.py`
- `scripts/lineage_prevalence.py`
- `scripts/mutation_profile.py`
- `scripts/reporting_lag.py`
- `scripts/resolve_lineage.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
