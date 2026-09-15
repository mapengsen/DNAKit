# Choosing the MACS2 Peak Threshold

ChIP-Atlas pre-computes peak calls at four MACS2 stringency levels. The
API exposes this as the `threshold` parameter, with the value equal to
the MACS2 `-10 × log10(p-value)` cutoff used when peaks were called.

## The four threshold values

| `threshold` | MACS2 score | Approx. p-value | Use case |
|-------------|-------------|-----------------|----------|
| **50** | -10·log₁₀(p) ≥ 50 | ~1e-5 | Default — balanced sensitivity/specificity |
| 100 | -10·log₁₀(p) ≥ 100 | ~1e-10 | High-confidence binding only |
| 200 | -10·log₁₀(p) ≥ 200 | ~1e-20 | Very high — strongest binding |
| 500 | -10·log₁₀(p) ≥ 500 | ~1e-50 | Extreme — broad histone marks, super-enhancers |

Higher threshold ⇒ fewer peaks ⇒ higher fold enrichment when something
is genuinely bound, but easier to miss real binding for moderately
expressed factors.

## How to choose

The default of `50` is right for the great majority of analyses. Raise
it deliberately:

- **Sanity-checking** — if a factor looks enriched at threshold 50, try
  100. If it still looks enriched, the call is robust. If it vanishes,
  the signal probably came from low-confidence peaks.
- **Histone marks with very broad domains** — H3K27me3, H3K36me3 — sit
  comfortably at threshold 100 or 200, where narrow noise peaks have
  been stripped out.
- **Avoiding well-studied-factor blanket coverage** — at threshold 500,
  only the strongest binding events for a TF survive, which can help
  separate near-saturated coverage of TP53/CTCF from genuinely strong
  enrichment in your gene set.

Lower than 50 is not supported by ChIP-Atlas.

## Effect on the contingency table

Each threshold corresponds to a different peak set, which changes the
numerators in the Fisher exact contingency table:

- `a` (input genes with peaks) shrinks at higher thresholds.
- `c` (background genes with peaks) shrinks faster (the background pool
  is much larger), so fold enrichment for genuine signals rises.
- p-values may go up or down depending on the trade-off — empirically,
  factors with strong, consistent binding remain significant at higher
  thresholds; factors that only ever showed marginal signal drop out.

## Quick robustness check

Re-run the same gene list at two thresholds and compare the top-N factor
tables:

```python
default = run_enrichment_workflow(gene_list=genes, threshold=50, ...)
strict  = run_enrichment_workflow(gene_list=genes, threshold=100, ...)
shared  = set(default["enrichment_results"]["antigen"].head(20)) & \
          set(strict["enrichment_results"]["antigen"].head(20))
```

If `len(shared) >= 15` the call is stable; if it drops below 10, treat
the top-N at threshold 50 as exploratory rather than a robust list.

## Citation

Zhang Y. *et al.* (2008). Model-based analysis of ChIP-seq (MACS).
*Genome Biology* 9(9):R137.
[doi:10.1186/gb-2008-9-9-r137](https://doi.org/10.1186/gb-2008-9-9-r137)
