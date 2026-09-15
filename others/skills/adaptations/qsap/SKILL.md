---
name: qsap
description: >
  QSAP — quorum-sensing (QS) related gene annotation pipeline against the QSP
  database (DIAMOND and/or HMMER). Use when mining autoinducer / QS circuitry
  (AHL, AIP, AI-2-related families as covered by QSP) from proteins or
  assemblies. Upstream: https://github.com/chunxiao-dcx/QSAP (QSP DB:
  https://github.com/chunxiao-dcx/QSP). Route via microbial-mining. Related
  secondary metabolites → antismash; Fe siderophores → fegenie.
license: MIT
category: analysis-tools
tags: [QSAP, quorum-sensing, QS, autoinducer, signaling, genome-mining]
upstream: https://github.com/chunxiao-dcx/QSAP
stage: mining
---

## DNAKit execution

Use `dnakit skills doctor qsap` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# QSAP (quorum sensing)

Upstream: [chunxiao-dcx/QSAP](https://github.com/chunxiao-dcx/QSAP) ·
QSP database: https://github.com/chunxiao-dcx/QSP

## Citation

Pipeline/DB: Dai, C. QSAP / QSP — quorum-sensing protein resources
(https://github.com/chunxiao-dcx/QSAP). Cite the upstream repository and QSP
release used; add peer-reviewed QSP/QSAP citation when available in Methods.

See also [docs/references.md](references/upstream-docs/references.md).

## Analytical thinking

**Quorum sensing (群体感应)** coordinates population behaviors via signal
synthases/receptors/transporters. QSAP classifies QS-related sequences from
the QSP collection. Treat hits as genetic potential — not measured autoinducer
concentrations.

| Related claim | Prefer |
|---------------|--------|
| QS genes | `qsap` |
| Specialized metabolite / AHL-adjacent BGCs | `antismash` |
| Siderophores (Fe acquisition signals/competition) | `fegenie` |
| Broad metabolism | `dram` / `element-cycling` |

## How to run

```bash
# install deps per QSAP README; example pattern
perl QSAP.pl -i input.list -o qsap_out -s sub -n 8
# strategies: diamond / hmmscan / sub|union — follow upstream help
```

## Decision tree

```text
Cell–cell signaling / QS?
├─ QS gene annotation → qsap
├─ BGC / natural products → antismash
├─ Siderophores → fegenie
└─ Plasmid-borne QS context → genomad / mob-suite ∩ qsap
```

## Related skills

`antismash` · `fegenie` · `dram` · `genomad` · `mob-suite` ·
`microbial-mining` · `tool-selection`

<!-- DNAKit documentation-link portability edit, 2026-09-12: bundled repository docs locally; see UPSTREAM.md. -->
