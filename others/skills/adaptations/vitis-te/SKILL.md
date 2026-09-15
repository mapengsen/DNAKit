---
name: vitis-te
description: >
  Transposable-element annotation and curation for plant genomes (Vitis).
  Use for EDTA, TEtrimmer, TEsorter, curated libraries, RepeatMasker, LAI.
  Do not hook a classifier-only tool as a genome-wide scanner.
---

## DNAKit execution

Use `dnakit skills doctor vitis-te` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# vitis-te

EDTA is the usual de novo starter. It is not a finished library.

## Order

1. EDTA per assembly (`--species others` unless it is rice/maize)
2. TEtrimmer for boundaries
3. TEsorter for lineage
4. Manual spot-check (LTR false positives, LINE/SINE misses, CDS contamination)
5. Curated lib: drop CDS → CD-HIT ~85% → 80-80 collapse
6. Re-annotate with `--curatedlib` and/or RepeatMasker
7. LTR age + LAI only after the lib is curated

## Red lines

- Do not treat EDTA raw output as a gold-standard TE library.
- Classifier-only tools (e.g. homology classifiers) are for **naming**, not for scanning a genome or replacing trim.
- Do not publish another group's raw TE calls as yours.

## EDTA starter

```bash
EDTA.pl --genome assembly.fa --species others --sensitive 1 --anno 1 --threads 32
```

Related: `vitis-pangenome`, `bioinfo-pitfalls`.
