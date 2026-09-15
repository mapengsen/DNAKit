# Target Genes — Server Layout and TSV Schema

This reference describes where the pre-computed target-gene files live,
how they are named, and how to read the wide-format TSV the downloader
parses.

## URL pattern

```
https://chip-atlas.dbcls.jp/data/{genome}/target/{Protein}.{Distance}.tsv
```

| Token | Allowed values | Notes |
|-------|---------------|-------|
| `{genome}` | `hg38`, `hg19`, `mm10`, `mm9`, `rn6`, `dm6`, `dm3`, `ce11`, `ce10`, `sacCer3` | Same set as the WABI API |
| `{Protein}` | Case-sensitive gene symbol (`TP53`, `CTCF`, `MYC`) | Histone marks return 404 |
| `{Distance}` | `1`, `5`, or `10` | ±kb from TSS |

Examples:

- `https://chip-atlas.dbcls.jp/data/hg38/target/TP53.5.tsv` (TP53 ±5 kb)
- `https://chip-atlas.dbcls.jp/data/mm10/target/Ctcf.10.tsv` (mouse CTCF, ±10 kb)
- `https://chip-atlas.dbcls.jp/data/hg38/target/E2F1.1.tsv` (E2F1 ±1 kb)

## HEAD check before download

```
HEAD https://chip-atlas.dbcls.jp/data/{genome}/target/{Protein}.{Distance}.tsv
```

- HTTP 200 → data exists; safe to GET.
- HTTP 404 → no data. Common reasons:
  - Symbol case (`tp53` vs `TP53`) — names are case-sensitive.
  - Histone marks — not served by the target-genes endpoint; use the
    enrichment workflow.
  - The factor has no ChIP-seq experiments in ChIP-Atlas.

`download_target_genes.check_antigen_available()` performs this check
and prints actionable diagnostics.

## TSV schema (wide format)

Column order:

| # | Header pattern | Type | Meaning |
|---|----------------|------|---------|
| 1 | `Target_genes` | str | Gene symbol (one row per gene) |
| 2 | `{Protein}|Average` | float | Mean MACS2 score across all experiments |
| 3..N-1 | `{SRX_ID}|{cell_type}` | float | Per-experiment MACS2 score |
| N | `STRING` | float | Pre-embedded STRING interaction score |

Example header (TP53 ±5 kb, hg38):

```
Target_genes\tTP53|Average\tSRX21888064|22Rv1\tSRX21888065|22Rv1\t...\tSTRING
```

Properties to know about:

- Rows are pre-sorted by the Average column, descending. The first row
  is the strongest "consensus" target.
- Scores are MACS2 `-10·log10(Q)` values (see
  `macs2_binding_scores.md`).
- When a row has multiple peaks for one gene, ChIP-Atlas keeps the
  highest score.
- Column count varies wildly: TP53 has ~395 columns; less-studied
  factors might have <10.
- Row count varies: TP53 ~16,000 genes; smaller TFs <1,000.
- Cell-type names follow ChIP-Atlas conventions — spelling may differ
  from RNA-seq sources (e.g. `K-562` not `K562`).

## Related files on the server

Each antigen has a small family of related files. The skill uses only
the first; the others are documented here for the curious:

- `{Protein}.{1,5,10}.tsv` — target gene data at different distances.
- `STRING_{Protein}.{1,5,10}.tsv` — STRING-enhanced predictions in a
  separate file; the regular TSV already carries STRING scores so this
  file is rarely needed.
- `{SRX_ID}.{1,5,10}.tsv` — single-experiment target lists. Useful if
  you want one specific dataset rather than the consensus.

## How target genes are computed

ChIP-Atlas overlaps experiment-level MACS2 peak intervals with the TSS
of every protein-coding gene from UCSC refFlat using `bedtools window`.
Gene-name normalization follows official nomenclature (HGNC for human,
MGI for mouse, etc.). When multiple peaks land in the same window, only
the highest-scoring peak contributes to that gene's score for that
experiment.

## Citation

Zou Z. *et al.* (2024). ChIP-Atlas 3.0: a data-mining suite to explore
chromosome architecture together with large-scale regulome data.
*Nucleic Acids Research* 52(D1):D45–D53.
[doi:10.1093/nar/gkad884](https://doi.org/10.1093/nar/gkad884)
