# DNAKit integration

Canonical skill: `causal-genomics-heritability-partitioning`.

```bash
python -m pip install dnakit
dnakit skills doctor causal-genomics-heritability-partitioning
dnakit skills install causal-genomics-heritability-partitioning --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts causal-genomics-heritability-partitioning
dnakit skills run causal-genomics-heritability-partitioning --script examples/ldak_sumher.sh -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `examples/ldak_sumher.sh`
- `examples/ldsc_partitioned_h2.sh`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
