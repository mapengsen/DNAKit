# DNAKit integration

Canonical skill: `alignment-msa-statistics`.

```bash
python -m pip install dnakit
dnakit skills doctor alignment-msa-statistics
dnakit skills install alignment-msa-statistics --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts alignment-msa-statistics
dnakit skills run alignment-msa-statistics --script examples/capra_singh_jsd.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `examples/capra_singh_jsd.py`
- `examples/conservation_profile.py`
- `examples/entropy_analysis.py`
- `examples/gap_statistics.py`
- `examples/identity_matrix.py`
- `examples/kimura_protein_distance.py`
- `examples/pssm.py`
- `examples/substitution_counts.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
