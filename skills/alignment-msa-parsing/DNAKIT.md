# DNAKit integration

Canonical skill: `alignment-msa-parsing`.

```bash
python -m pip install dnakit
dnakit skills doctor alignment-msa-parsing
dnakit skills install alignment-msa-parsing --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.


## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts alignment-msa-parsing
dnakit skills run alignment-msa-parsing --script examples/a2m_a3m_io.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `examples/a2m_a3m_io.py`
- `examples/analyze_alignment.py`
- `examples/clean_alignment.py`
- `examples/consensus_sequence.py`
- `examples/find_conserved.py`
- `examples/gap_analysis.py`
- `examples/henikoff_weights.py`
- `examples/mi_apc.py`
- `examples/neff.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
