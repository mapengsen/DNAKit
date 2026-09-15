---
name: vitis-popgen
description: >
  Population-genetics workflow from a clean VCF: clones, PCA, ADMIXTURE,
  Fst, TreeMix, Dsuite, selection. Use for Vitis / perennial crops where
  clones and families wreck structure plots.
---

## DNAKit execution

Use `dnakit skills doctor vitis-popgen` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# vitis-popgen

Grapevine is clonal. Filter relatives before you interpret clusters.

## Order

```
VCF (already QC'd)
  → LD prune          plink --indep-pairwise
  → kinship / clones  KING or PLINK --genome
  → drop clones / parent-child / full sibs
  → PCA               GCTA / PLINK
  → ADMIXTURE K=2..   + align across K (pong)
  → optional tree     IQ-TREE / FastTree
  → Fst / π / ROH
  → gene flow         TreeMix, Dsuite, ADMIXTOOLS2
  → selection         only after structure is honest
```

## Red lines

1. No HWE / missingness / clone filter → do not publish ADMIXTURE.
2. Do not put highly related samples into the learning set for ADMIXTURE/PCA.
3. Project aDNA or sparse chip samples; do not let them define K.
4. Multiple testing on genome-wide scans is not optional.

## Clonal crops

- One accession, many cuttings → treat as one individual.
- Secondary centres / feral vines: label geography and use, not just country.

Related: `grapevine-adna`, `grapeancestry`, `bioinfo-pitfalls`.
