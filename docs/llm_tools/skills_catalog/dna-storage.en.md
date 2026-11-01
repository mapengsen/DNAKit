# DNA data storage

**1** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`dnakit-dna-storage`](#dnakit-dna-storage) | dnakit | Encode bytes into DNA payloads and decode intact packets through the optional Chamaeleo Goldman codec. Use for computational DNA storage demonstrations and cod… |

## `dnakit-dna-storage` {#dnakit-dna-storage}

Encode bytes into DNA payloads and decode intact packets through the optional Chamaeleo Goldman codec. Use for computational DNA storage demonstrations and codec integration.

- **Source:** [DNAKit contributors; method authors cited in SKILL.md and DNAKIT.md @ `dnakit-`](https://github.com/mapengsen/DNAKit/tree/main/skills/dnakit-dna-storage)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `Chamaeleo`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info dnakit-dna-storage
dnakit skills install dnakit-dna-storage --agent claude-code   # or --agent codex / cursor
```
