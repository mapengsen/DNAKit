# MACS2 Binding Scores

The scores in ChIP-Atlas target-gene TSVs come from MACS2 peak calls
and follow a single, simple transformation. This reference explains the
score scale and why the *averaged* score in the ranking columns has a
different interpretation from a single experiment's score.

## Score formula

For every peak call, MACS2 emits a Benjamini-Hochberg-corrected
q-value. ChIP-Atlas reports:

$$\text{score} = -10 \times \log_{10}(\text{MACS2 Q-value})$$

Higher score = stronger binding evidence. The mapping:

| Score | MACS2 Q-value | What it means in a single experiment |
|-------|--------------|--------------------------------------|
| 0 | no peak | No binding detected |
| 50 | 1e-5 | Default MACS2 threshold (weak but detectable) |
| 100 | 1e-10 | Moderate binding |
| 200 | 1e-20 | Strong binding |
| 500 | 1e-50 | Very strong binding |
| 1000 | 1e-100 | Extremely strong binding |

## Handling multiple peaks per gene

When more than one peak from a single experiment overlaps a gene's TSS
window, ChIP-Atlas keeps only the **highest** score for that pair. This
means each experiment contributes at most one score per gene.

## The `{Protein}|Average` column

The Average column is the arithmetic mean of MACS2 scores across all
experiments for that factor:

$$\text{average} = \frac{1}{N_{\text{experiments}}} \sum_i \text{score}_i$$

Useful interpretive bins:

| Average | Reading |
|---------|---------|
| > 500 | Consistently strong binding across many experiments — high-confidence target |
| 100–500 | Strong, but possibly cell-type-specific |
| 50–100 | Moderate; check binding rate before believing |
| < 50 | Weak; likely cell-type-specific or noisy |

## Why the average dilutes

Because the average runs across **all** experiments — including ones
where the factor was *not* bound to this gene — every zero pulls the
mean down. A gene with score 1000 in 5 of 395 experiments has an
average of `(5 × 1000 + 390 × 0) / 395 ≈ 12.7`, not 1000.

This has two practical consequences:

1. The score-to-Q-value table above does **not** apply to the average
   column. An average of 500 means the *consensus* reaches that level,
   not that every experiment hit Q ≤ 1e-50.
2. Cell-type-specific targets get lower averages than universal targets,
   even when their binding is very strong in the relevant cell type.

## Complementary metrics

The skill computes two additional summary columns to recover the
information the average dilutes away:

| Column | Definition |
|--------|------------|
| `max_score` | Highest score in any single experiment. Captures strong cell-type-specific binding the average obscures. |
| `binding_rate` | Fraction of experiments (0–1) with any non-zero score. Higher rate = more consistent across cell types. |
| `num_bound` | Raw count of experiments with non-zero score. Pair this with `num_experiments` (the denominator) for a quick `n/N`. |

A useful pattern: high `max_score` + low `binding_rate` flags a strong
but cell-type-specific target; high `max_score` + high `binding_rate`
flags a universal one.

## Re-subsetting by cell type

When you pass `cell_types=["MCF-7"]` (or any list) to the workflow, the
filter recomputes `avg_score`, `max_score`, `binding_rate`, and
`num_bound` using **only** the matching experiments. The wide-format
`experiment_data` is also subset to those columns. This is the cleanest
way to get a cell-type-specific ranking.

## Peak calling parameters

ChIP-Atlas runs `MACS2 callpeak` with default settings:

- Q-value threshold: `1e-5` (any weaker peaks never enter the data)
- Genome size: `-g hs` (human), `-g mm` (mouse), etc.
- Tag size: inferred from the FASTQ

Only peaks passing the `1e-5` threshold appear in the target-gene
tables. Going below score 50 in the ranking means the gene was bound
just above this floor in a few experiments and zero everywhere else.

## Citations

Zhang Y. *et al.* (2008). Model-based analysis of ChIP-seq (MACS).
*Genome Biology* 9(9):R137.
[doi:10.1186/gb-2008-9-9-r137](https://doi.org/10.1186/gb-2008-9-9-r137)

Zou Z. *et al.* (2024). ChIP-Atlas 3.0: a data-mining suite to explore
chromosome architecture together with large-scale regulome data.
*Nucleic Acids Research* 52(D1):D45–D53.
[doi:10.1093/nar/gkad884](https://doi.org/10.1093/nar/gkad884)
