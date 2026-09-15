# DNAKit integration

Canonical skill: `recombinator`.

```bash
python -m pip install dnakit
dnakit skills doctor recombinator
dnakit skills install recombinator --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts recombinator
dnakit skills run recombinator --script recombinator.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `recombinator.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.

## Simulation data

DNAKit bundles the attributed inheritance engine. Use `--demo` for a synthetic one-locus example. For supplied parent genomes, use `--data-dir` with `GENOMES/`, `trait_registry.json` and `disease_registry.json`, or explicit registry paths. The original historical-person demo genomes were absent from the pinned repository and are not fabricated here. The Python API is `dnakit.workflows.simulate_inheritance(...)`, with optional offspring count, mutation rate and seed. This is independent allele inheritance, not chromosome-scale crossover or a clinical predictor.
