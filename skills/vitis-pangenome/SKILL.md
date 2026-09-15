---
name: vitis-pangenome
description: >
  Vitis / crop pangenome decisions: when to build a graph, PAV, reference
  choice for mixed cohorts and chips. Use for minigraph-cactus vs PGGB,
  PanGenie, or "should we build another genus graph".
---

## DNAKit execution

Use `dnakit skills doctor vitis-pangenome` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# vitis-pangenome

A new genus-wide graph is rarely the next paper. Most jobs are reference choice, PAV, and genotyping.

## Decide first

1. **Question** — presence/absence, SV-eQTL, or just mapping a chip cohort?
2. **Already-published graphs** — reuse and genotype; do not rebuild for prestige.
3. **Mixed modern + wild + aDNA** — a linear reference + projection often beats a half-built graph.

## Tool defaults

| Job | First tool |
|-----|------------|
| Reference-guided graph | minigraph-cactus |
| Unbiased all-vs-all | PGGB (expensive) |
| Genotype a cohort on a graph | vg / PanGenie |
| Pairwise SV / synteny | SyRI |
| Chip / mixed-cohort mapping | pick one linear ref (document why), not a moving graph |

## Red lines

- Do not invent a "super-pangenome" from a handful of assemblies.
- Do not mix graph coordinates and chip coordinates without a lift.
- Field move in crops is mixed-variant / dosage GWAS, not another genus GFA.

Related: `grapevine-chip`, `vitis-popgen`, `grapevine-adna`.
