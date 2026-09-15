# DNAKit integration

Canonical skill: `reverse-complement`.

```bash
python -m pip install dnakit
dnakit skills doctor reverse-complement
dnakit skills install reverse-complement --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts reverse-complement
dnakit skills run reverse-complement --script examples/basic_operations.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `examples/basic_operations.py`
- `examples/iupac_and_features.py`
- `examples/palindrome_check.py`
- `examples/primer_design.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
