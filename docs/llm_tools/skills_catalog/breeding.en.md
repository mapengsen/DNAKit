# Genomic breeding and agriculture

**1** Skills. A Skill is a set of instructions that guides an Agent through an analysis; the software, reference data and APIs it uses are separate dependencies. See [Skills](../skills.md) for searching, installing and using them; the Skills tab of the web demo shows each Skill in full.

| Skill | Source | Summary |
| --- | --- | --- |
| [`dnakit-genomic-breeding`](#dnakit-genomic-breeding) | dnakit | Compute diploid genomic relationships and fixed-lambda RR-BLUP predictions with DNAKit. Use for genotype-to-phenotype breeding research when aligned dosage mat… |

## `dnakit-genomic-breeding` {#dnakit-genomic-breeding}

Compute diploid genomic relationships and fixed-lambda RR-BLUP predictions with DNAKit. Use for genotype-to-phenotype breeding research when aligned dosage matrices and training phenotypes are available.

- **Source:** [DNAKit contributors; method authors cited in SKILL.md and DNAKIT.md @ `dnakit-`](https://github.com/mapengsen/DNAKit/tree/main/skills/dnakit-genomic-breeding)
- **License:** `MIT`
- **Runs with:** Python modules `dnakit`, `numpy`
- **Checks:** metadata passed, resources checksummed, execution not_run

```bash
dnakit skills info dnakit-genomic-breeding
dnakit skills install dnakit-genomic-breeding --agent claude-code   # or --agent codex / cursor
```
