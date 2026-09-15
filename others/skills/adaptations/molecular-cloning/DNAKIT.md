# DNAKit integration

Canonical skill: `molecular-cloning`.

```bash
python -m pip install dnakit
dnakit skills doctor molecular-cloning
dnakit skills install molecular-cloning --agent codex
```

Follow SKILL.md for the scientific procedure and interpretation. DNAKit supplies the skill and its listed local resources. Real scientific libraries, command-line programs, models, reference databases and service credentials retain their own requirements.

## Database functions

```python
import dnakit
# Inspect the exact parameters before calling a service:
schema = dnakit.databases.describe('dna_gibson_design')
```

DNAKit bundles the attributed adapter implementations; installing ToolUniverse is unnecessary. Provider result fields are preserved. DNAKit raises exceptions on provider errors; handle these without inventing missing evidence.

| Upstream operation | DNAKit function |
| --- | --- |
| `DNA_gibson_design` | `dnakit.databases.dna_gibson_design` |
| `DNA_golden_gate_assemble` | `dnakit.databases.dna_golden_gate_assemble` |
| `DNA_golden_gate_design` | `dnakit.databases.dna_golden_gate_design` |
| `DNA_virtual_digest` | `dnakit.databases.dna_virtual_digest` |

## Bundled scripts

Find scripts and run the relevant one with its documented arguments:

```bash
dnakit skills scripts molecular-cloning
dnakit skills run molecular-cloning --script scripts/cloning_qc.py -- <documented arguments>
```

The command resolves files from the installed DNAKit package; a checkout of the upstream skill collection is unnecessary. It forwards standard input/output and exit status. Python callers can use `dnakit.workflows.run(skill_id, script, arguments=[...], workdir=...)` for bounded captured output.

- `scripts/cloning_qc.py`

## Provenance and validation

Original authors, license and revision remain in [UPSTREAM.md](UPSTREAM.md). DNAKit adaptation changes are recorded in the catalog. Packaging and interface checks do not establish scientific validity or completion of the full workflow.
