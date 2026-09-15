# DNAKit integration

Canonical skill: `methylation-aggregation`.

```bash
python -m pip install dnakit
dnakit skills doctor methylation-aggregation
dnakit skills install methylation-aggregation --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts methylation-aggregation
dnakit skills run methylation-aggregation --script scripts/validate_methylation.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `scripts/validate_methylation.py`

## ENCODE

Use `dnakit.encode` for ENCODE queries, downloads and provenance tracking. The API retains the original options; `encode_search_experiments(...)` is now `dnakit.encode.search_experiments(...)`. It returns the original JSON text. Credential storage uses `dnakit[encode-credentials]`; public searches need no credentials.


## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
