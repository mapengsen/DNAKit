# Promoter / TSS Window Mechanics

This file explains how the API decides which peaks count as "near" each
gene, and what the optional Ensembl coordinate verification step does.

## The API does the gene-to-region mapping for you

When you submit with `typeA=gene`, the ChIP-Atlas server walks each
symbol through its own internal RefSeq annotation. The window applied is
controlled by:

- `distanceUp` — bases upstream of the TSS (default 5000)
- `distanceDown` — bases downstream of the TSS (default 5000)

The default 10 kb window is broad — fine for asking "any nearby
regulator" but loose enough to capture distal enhancer hits. Narrower
windows tighten the test toward the core promoter.

| Window preset | distanceUp | distanceDown | Total | Best for |
|---------------|-----------|--------------|-------|----------|
| Default (API) | 5000 | 5000 | 10 kb | Broad regulatory scan |
| Standard | 2000 | 500 | 2.5 kb | Core promoter + immediate flank |
| Narrow | 500 | 100 | 600 bp | Core promoter only |
| Wide | 10000 | 10000 | 20 kb | Distal enhancers / open chromatin |

Pick a preset based on what you mean by "promoter"; don't fight the API
default for arbitrary reasons.

## TSS direction depends on strand

If you ever need to compute regions yourself (e.g. for an offline
overlap), follow these conventions:

| Strand | TSS coordinate | Region |
|--------|---------------|--------|
| `+` | `gene_start` | `[TSS - distanceUp, TSS + distanceDown]` |
| `−` | `gene_end` | `[TSS - distanceDown, TSS + distanceUp]` |

Symbolically: the upstream window is always 5′ of the TSS, regardless of
strand.

## Ensembl coordinate verification (optional)

`run_enrichment_workflow.py` calls `convert_genes_to_regions.py` after
the API returns. This populates the `input_regions` list used by
downstream visualization scripts. It is *not* used to re-run the
enrichment — those numbers come from ChIP-Atlas's own RefSeq mapping.

The Ensembl lookup uses:

```
POST https://rest.ensembl.org/lookup/symbol/{species}
Body: {"symbols": ["TP53", "CDKN1A", "BAX"]}
```

with retry-and-fallback to per-symbol GETs.

Sequential per-gene retries are throttled to 0.5 seconds and back off
exponentially on HTTP 429.

### Common Ensembl pitfalls

- **Different gene database than ChIP-Atlas.** The API uses RefSeq;
  Ensembl uses Ensembl. Their symbol coverage overlaps but isn't
  identical — a gene can be resolved by ChIP-Atlas and fail in Ensembl
  (or vice versa).
- **Retired HGNC symbols.** The converter applies a small built-in alias
  map (IL8→CXCL8, TNFSF6→FASLG, MLL→KMT2A, etc.). Other retired symbols
  may simply fail. Surface the failures in the summary so users know to
  update.
- **Transient outages.** When Ensembl is slow or down, the workflow
  retries once after a 60-second wait. If that fails, the enrichment is
  reported anyway with `input_regions=[]` and a prominent caveat in the
  markdown summary.

### What discrepancy means

The summary report distinguishes three cases:

| Ensembl mapped | API analysed | Reading |
|----------------|--------------|---------|
| 0 / N | N | Ensembl is unavailable but enrichment is valid. Loss is "no independent coordinate verification." |
| < N / N | N | Mix of resolution issues. Cite both counts. |
| N / N | N | Both systems agree. No caveat needed beyond data-availability bias. |

Do not dismiss a 0/N Ensembl outcome as "optional." It just means the
two systems disagree about gene-name resolution and the enrichment
relies entirely on the API's RefSeq mapping.

## Citation

Yates A.D. *et al.* (2020). Ensembl 2020. *Nucleic Acids Research*
48(D1):D682–D688.
[doi:10.1093/nar/gkz966](https://doi.org/10.1093/nar/gkz966)

RefSeq: Pruitt K.D. *et al.* (2014). RefSeq: an update on mammalian
reference sequences. *Nucleic Acids Research* 42(D1):D756–D763.
[doi:10.1093/nar/gkt1114](https://doi.org/10.1093/nar/gkt1114)
