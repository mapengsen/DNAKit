# STRING Interaction Scores

The `STRING` column in the ChIP-Atlas target-gene TSV reflects evidence
from the STRING database. STRING is queried once during the ChIP-Atlas
data pipeline; the column is pre-embedded, so the skill does NOT make a
separate STRING API call.

## How the scores enter the TSV

For each (protein, gene) pair, ChIP-Atlas keeps the STRING score when
ALL of these hold:

- `item_id_a` matches the queried antigen (the TF / protein).
- `item_id_b` matches the target gene.
- Interaction mode is `"expression"` (regulatory direction).
- `a_is_acting = 1` (the antigen is acting on the gene).

If none of the recorded interactions meet these criteria, the column
value is `0`.

## Score interpretation

STRING scores run 0–1000 in standard STRING terms. The ChIP-Atlas
column preserves that scale. Reference bands:

| Score | STRING confidence band |
|-------|-----------------------|
| 0 | No supporting evidence in STRING |
| 1–399 | Low confidence |
| 400–699 | Medium confidence |
| 700–899 | High confidence |
| 900–1000 | Highest confidence |

ChIP-Atlas applies a combined-score threshold internally; very low
STRING scores are unlikely to appear at all.

## Combining STRING with binding

The strongest predictions have **both**: high MACS2 binding evidence
*and* high STRING evidence. A useful rubric:

| Binding (avg) | STRING | Confidence |
|---------------|--------|------------|
| High (>200) | High (>400) | Very high — direct evidence + regulatory match |
| High (>200) | 0 | High — direct binding; novel target or STRING gap |
| Low (<50) | High (>400) | Moderate — known interaction, weak / indirect binding |
| Low (<50) | 0 | Low — minimal evidence |

## "STRING = 0" ≠ "not a target"

This is the single most important caveat. STRING coverage is sparse and
strongly biased toward well-studied interactions. A score of 0 means
STRING lacks evidence for the specific (protein, gene) pair — not that
the pair is biologically irrelevant.

Concrete example: BBC3 (PUMA) is one of the best-characterised TP53
apoptosis effectors. In the ChIP-Atlas TP53 target-gene file, its
STRING score is 0, even though its MACS2 binding is strong. The
0 reflects a coverage gap, not biology.

When summarizing results for a user, always note: **absence of STRING
evidence is not evidence of absence**. Many genuine targets have STRING
score 0.

## Limitations to call out

- STRING "expression" mode aggregates co-expression, text-mining,
  experimental, and database evidence. A high score can reflect
  protein-protein interaction or pathway co-membership rather than
  direct transcriptional regulation.
- STRING's combined-score threshold is conservative; novel
  cell-type-specific regulatory relationships are unlikely to be
  represented even when they exist in the literature.
- The TSV column is a snapshot from when ChIP-Atlas built the data;
  STRING evolves between releases.

## Citations

Szklarczyk D. *et al.* (2023). The STRING database in 2023.
*Nucleic Acids Research* 51(D1):D638–D646.
[doi:10.1093/nar/gkac1000](https://doi.org/10.1093/nar/gkac1000)

Zou Z. *et al.* (2024). ChIP-Atlas 3.0: a data-mining suite to explore
chromosome architecture together with large-scale regulome data.
*Nucleic Acids Research* 52(D1):D45–D53.
[doi:10.1093/nar/gkad884](https://doi.org/10.1093/nar/gkad884)
