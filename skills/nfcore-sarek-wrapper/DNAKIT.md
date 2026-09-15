# DNAKit integration

Canonical skill: `nfcore-sarek-wrapper`.

```bash
python -m pip install dnakit
dnakit skills doctor nfcore-sarek-wrapper
dnakit skills install nfcore-sarek-wrapper --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts nfcore-sarek-wrapper
dnakit skills run nfcore-sarek-wrapper --script _isolated_imports.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `_isolated_imports.py`
- `command_builder.py`
- `errors.py`
- `executor.py`
- `nfcore_sarek_wrapper.py`
- `outputs_parser.py`
- `params_builder.py`
- `pipeline_source.py`
- `preflight.py`
- `provenance.py`
- `remap_paths.py`
- `reporting.py`
- `samplesheet_builder.py`
- `schemas.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
