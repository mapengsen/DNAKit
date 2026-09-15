# Statistical Foundation of the Enrichment Call

This reference explains the math behind the per-experiment q-values and
fold enrichments so you can interpret edge cases correctly.

## Fisher's exact test

For every experiment in the database, the server fills a 2×2 table:

|                         | Peak overlap | No peak overlap |
|-------------------------|:------------:|:---------------:|
| **Input gene set**      | a            | b               |
| **RefSeq background**   | c            | d               |

where

- `a` — input genes with at least one MACS2 peak in this experiment's
  TSS window
- `b` — input genes without peaks
- `c` — RefSeq genes with peaks (excluding the input set)
- `d` — RefSeq genes without peaks
- `a + b = n_input`, `c + d ≈ 18,600` for human (varies by genome)

The API reports a **one-sided** test for enrichment (odds ratio > 1) and
emits `log10(p-value)` to preserve precision near zero.

## Fold enrichment

$$\text{FE} = \frac{a / n_{input}}{c / n_{background}}$$

Useful interpretive bins:

| Fold enrichment | Interpretation |
|-----------------|----------------|
| > 10× | Very strong — likely a direct regulator |
| 5–10× | Strong |
| 2–5× | Moderate |
| < 2× | Weak; may be background |
| < 1× | Depletion (API reports ~1e-6) |
| > 100,000× | **Sentinel** — see below |

## Sentinel values

When the RefSeq background overlap count `c` is zero, the FE formula
divides by zero. The API emits a sentinel of roughly `1,000,000`.

How to recognise a sentinel:

- FE ≥ 100,000×
- Small `regions_with_overlaps` (1 / N is typical)
- Often `q_value > 0.05` (the Fisher exact still penalises for `a` being
  small)

How the skill handles them:

- Ranking is by q-value; sentinels don't take over the top because they
  rarely reach significance.
- `enrichment_results_top20.csv` and the visualization panels require
  `regions_with_overlaps >= 2` to drop these.
- Volcano / scatter plots cap FE at 1000× for display so the axis is
  readable.

## Benjamini-Hochberg correction

Each batch of results undergoes BH correction across all tested
experiments. The procedure:

1. Sort by ascending p-value.
2. Compute `q_i = p_i * N / rank_i`.
3. Enforce monotonicity from the bottom up.

Typical thresholds:

- `q < 0.001` — highly significant after correction
- `q < 0.01` — significant
- `q < 0.05` — genome-wide significant
- `q ≥ 0.05` — not significant

The API reports `log10(q-value)`; the parser converts back to raw
probabilities for downstream calculations.

## Overlap rate

A simple ratio, reported in the `overlap_rate` column:

$$\text{overlap rate} = \frac{\text{regions\_with\_overlaps}}{\text{total\_regions}}$$

Helpful framing:

- > 50% — core regulon
- 20–50% — key targets
- < 20% — subset or indirect

Cite the median over significant experiments rather than picking a single
extreme value, and prefer the aggregated factor table for medians across
all of a factor's experiments.

## Data-availability bias

ChIP-Atlas mirrors the public ChIP-seq universe, which is far from
uniform. Two structural biases that *always* deserve mention:

1. **Well-studied TFs are over-represented.** TP53, RELA, CTCF, MYC, and
   ESR1 each have hundreds to thousands of experiments. A factor's
   prominence in the top-K can partly reflect being well-studied.
   Surface the "Experiments" column from the aggregated factor table
   so the reader can judge.
2. **Common cell lines dominate.** K-562, HeLa, MCF-7, HepG2, GM12878,
   and a handful of others contribute disproportionately. Enrichment in
   those lines is more likely to be detected for purely statistical
   reasons.

Recommendations:

- Rank by q-value, not fold enrichment.
- For each top factor, compare `Experiments` vs `Sig` (the significant
  count). A factor with 137 experiments but only 5 significant is
  showing cell-type-specific enrichment, not universal regulation.
- For "consensus regulator" claims, validate with orthogonal evidence
  (expression correlation, perturbation, motif analysis).

## When the API returns fewer regions than genes submitted

Look at `api_total_regions` in the workflow output. If it is below
`len(input_genes)`, one or more genes were dropped. Likely causes:

- The gene is not in RefSeq (common for very recently renamed symbols).
- The gene has a duplicate / collapsed promoter region.
- The symbol failed to resolve to a coordinate.

The API does **not** report which specific genes were dropped. Surface
the discrepancy to the user, but do not speculate about which gene was
lost.

## Citation

Zou Z. *et al.* (2024). ChIP-Atlas 3.0: a data-mining suite to explore
chromosome architecture together with large-scale regulome data.
*Nucleic Acids Research* 52(D1):D45–D53.
[doi:10.1093/nar/gkad884](https://doi.org/10.1093/nar/gkad884)

Oki S. *et al.* (2018). ChIP-Atlas: a data-mining suite powered by full
integration of public ChIP-seq data. *EMBO Reports* 19(12):e46255.
[doi:10.15252/embr.201846255](https://doi.org/10.15252/embr.201846255)
