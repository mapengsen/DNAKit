# DNA nanotechnology

**1** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`dnakit-dna-nanotechnology`](#dnakit-dna-nanotechnology) | dnakit | Construct an assigned scadnano DNA duplex and export caller-supplied scadnano geometry into oxDNA topology and configuration text using DNAKit. |

## `dnakit-dna-nanotechnology` {#dnakit-dna-nanotechnology}

Construct an assigned scadnano DNA duplex and export caller-supplied scadnano geometry into oxDNA topology and configuration text using DNAKit.

- **Source:** [DNAKit contributors; method authors cited in SKILL.md and DNAKIT.md @ `dnakit-`](https://github.com/mapengsen/DNAKit/tree/main/skills/dnakit-dna-nanotechnology)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `scadnano`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info dnakit-dna-nanotechnology
dnakit skills install dnakit-dna-nanotechnology --agent claude-code   # or --agent codex / cursor
```
