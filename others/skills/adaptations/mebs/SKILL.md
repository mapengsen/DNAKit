---
name: mebs
description: >
  MEBS (Multigenomic Entropy-Based Score) — score genomes/metagenomes for
  biogeochemical cycle machinery via informative Pfam domains; strongest for
  sulfur (S), also supports N, O, CH4, and Fe cycle scores. Use for comparative
  cycle enrichment across large (meta)genome collections. Upstream:
  https://github.com/valdeanda/mebs (also eead-csic-compbio/metagenome_Pfam_score).
  Route via element-cycling. Gene-level Fe → fegenie; multi-element diagrams →
  metabolic.
license: MIT
category: analysis-tools
tags: [MEBS, sulfur, S, N, CH4, Fe, biogeochemistry, element-cycling]
upstream: https://github.com/valdeanda/mebs
stage: mining
---

## DNAKit execution

Use `dnakit skills doctor mebs` to inspect scientific dependencies. Read [DNAKIT.md](DNAKIT.md) for this skill's DNAKit functions and bundled script entrypoints. Keep the workflow, input checks and reporting requirements below.


# MEBS

Upstream: [valdeanda/mebs](https://github.com/valdeanda/mebs) ·
legacy/core: https://github.com/eead-csic-compbio/metagenome_Pfam_score

## Citation

De Anda, V. *et al.* MEBS, a software platform to evaluate large (meta)genomic
collections according to their metabolic machinery: unraveling the sulfur cycle.
*Gigascience* **6**, 1–17 (2017).
https://doi.org/10.1093/gigascience/gix096

See also [docs/references.md](references/upstream-docs/references.md).

## Analytical thinking

MEBS compresses curated cycle machinery into **relative-entropy scores** over
Pfam domains. Benchmarked deeply on the **sulfur (S)** cycle; documentation
also supports **N**, **O**, **CH₄**, and **Fe** scoring modes. Ideal for
ranking many genomes/metagenomes by cycle “machinery presence,” not for
drawing full METABOLIC-style pathway cartoons.

| vs | Prefer MEBS when |
|----|------------------|
| `metabolic` | Comparative S/N/CH₄/Fe scores across large collections |
| `fegenie` | Score-level Fe vs explicit Fe gene neighborhoods |
| `dram` | Cycle enrichment metric vs distilled product tables |

Input: directory of `.faa` peptide FASTA; use genomic vs metagenomic mode per
docs (fragmentary metagenomes need the metagenomic entropy tables).

## How to run

```bash
# clone / install per versioned README (v1.1 / v1.2 manuals differ)
# typical: mebs.pl -i pep_dir -t genomic|metagenomic …  (verify installed flags)
```

## Decision tree

```text
Cycle scoring need?
├─ S (or N/O/CH₄/Fe) entropy scores → mebs
├─ C/N/S diagrams + networks → metabolic
├─ Fe gene calls → fegenie
└─ Broad distillate → dram
```

## Related skills

`element-cycling` · `metabolic` · `fegenie` · `dram` · `microbial-mining` ·
`tool-selection`

<!-- DNAKit documentation-link portability edit, 2026-09-12: bundled repository docs locally; see UPSTREAM.md. -->
