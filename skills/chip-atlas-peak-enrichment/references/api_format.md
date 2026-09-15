# ChIP-Atlas Enrichment Analysis — API Specification

This is the request/response contract for the WABI Enrichment Analysis
endpoint. The scripts in this skill call it directly; this reference is
the place to look when something fails or you need to understand a column
in the raw response.

## Endpoint

```
POST https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/
```

The same host serves the Diff Analysis endpoint. The `antigenClass`
value selects mode: pass an antigen class (e.g. `"TFs and others"`) for
enrichment, or `"diffbind"`/`"dmr"` for diff analysis.

The job model is asynchronous: submit → poll → retrieve. Submits return
a `requestId`; subsequent GETs use that ID with `?info=status` and then
`?info=result&format=tsv`.

## POST parameters

All fields are required form data. Empty placeholders cannot be omitted.

| Field | Required value or example | Notes |
|-------|--------------------------|-------|
| `address` | `""` | Reserved field; must be present but empty |
| `format` | `"json"` | Enrichment responses are JSON (diff uses text) |
| `result` | `"www"` | Reserved field |
| `genome` | `hg38` / `hg19` / `mm10` / `mm9` / `rn6` / `dm6` / `dm3` / `ce11` / `ce10` / `sacCer3` | Ten supported assemblies |
| `antigenClass` | e.g. `"TFs and others"`, `"Histone"` | Must be non-empty |
| `cellClass` | `"All cell types"` for the unrestricted query, or a category name | Must be non-empty |
| `threshold` | `"50"` (default), `"100"`, `"200"`, `"500"` | MACS2 -10 log10(p) threshold |
| `typeA` | `"gene"` | Tells the API the payload is gene symbols |
| `bedAFile` | Newline-separated gene symbols | Sanitize to alphanumerics + underscore before submitting |
| `typeB` | `"refseq"` | RefSeq promoters as background |
| `bedBFile` | `"empty"` | Reserved placeholder for the background payload |
| `permTime` | `"1"` | Reserved |
| `title` | Free text | Alphanumerics + underscore — appears in logs |
| `descriptionA` | Free text | Display label for the input set |
| `descriptionB` | Free text | Display label for the background |
| `distanceUp` | `"5000"` (default) | Bases upstream of TSS |
| `distanceDown` | `"5000"` (default) | Bases downstream of TSS |

## Supported assemblies

| Organism | Assemblies |
|----------|-----------|
| *Homo sapiens* | hg38, hg19 |
| *Mus musculus* | mm10, mm9 |
| *Rattus norvegicus* | rn6 |
| *Drosophila melanogaster* | dm6, dm3 |
| *Caenorhabditis elegans* | ce11, ce10 |
| *Saccharomyces cerevisiae* | sacCer3 |

## Antigen-class vocabulary

| Value | Meaning |
|-------|---------|
| `"TFs and others"` | Transcription factors and co-factors (most common) |
| `"Histone"` | Histone post-translational modifications (H3K27ac, H3K4me3, etc.) |
| `"RNA polymerase"` | Pol II / Pol III variants |
| `"ATAC-Seq"` | Open-chromatin assays |
| `"DNase-seq"` | DNase hypersensitivity |
| `"Bisulfite-Seq"` | DNA methylation |
| `"Input control"` | Control experiments |
| `"Annotation tracks"` | Genomic annotations |

## Cell-class vocabulary (hg38)

`"All cell types"`, `"Adipocyte"`, `"Blood"`, `"Bone"`, `"Breast"`,
`"Cardiovascular"`, `"Digestive tract"`, `"Epidermis"`, `"Gonad"`,
`"Kidney"`, `"Liver"`, `"Lung"`, `"Muscle"`, `"Neural"`, `"Others"`,
`"Pancreas"`, `"Placenta"`, `"Pluripotent stem cell"`, `"Prostate"`,
`"Uterus"`.

Other genomes use the same names with smaller subsets — invalid values
return HTTP 400.

## Submit response

JSON body:

```json
{
  "requestId": "wabi_chipatlas_2025-0517-1942-12-987-654321",
  "current-time": "2025-05-17 19:42:12"
}
```

If `requestId` is missing, the body contains an `error-message` /
`Message` field with the failure reason.

## Polling

```
GET https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/{requestId}?info=status
```

Body is plain text with one informative line:

```
status: running
```

or

```
status: finished
```

Status `error`/`failed` indicates a server-side failure — read the
remainder of the body for diagnostic detail.

Enrichment jobs are typically 1–3 minutes. "All cell types" with hundreds
of experiments per cell class can stretch to 10 minutes; raise the poll
timeout accordingly.

## Result retrieval

```
GET https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/{requestId}?info=result&format=tsv
```

Response body is a headerless TSV with 11 columns, one row per
experiment.

| # | Column | Type | Example |
|---|--------|------|---------|
| 1 | `experiment_id` | str | `SRX21147142` |
| 2 | `antigen_class` | str | `TFs and others` |
| 3 | `antigen` | str | `TP53` |
| 4 | `cell_type_class` | str | `Blood` |
| 5 | `cell_type` | str | `MOLM-13` |
| 6 | `num_peaks` | int | `99` |
| 7 | `overlap_a` | str | `"4/5"` (input genes with peaks / input gene count) |
| 8 | `overlap_b` | str | `"10/18851"` (background hits / background size) |
| 9 | `log10_pvalue` | float | `-12.0222` |
| 10 | `log10_qvalue` | float | `-8.15177` |
| 11 | `fold_enrichment` | float | `1508.08` |

Caveats:

- The API encodes very-large values as a sentinel near `1e6` when the
  background overlap is zero (see `enrichment_statistics.md`).
- `log10_pvalue` / `log10_qvalue` of `0` means `p=1` (non-significant);
  the parser converts to raw probabilities.
- Sort order: rows are typically returned by ascending q-value.

## Common failure modes

- **400 Bad Request** — almost always a missing or unrecognised
  `antigenClass` / `cellClass`. Both must be non-empty.
- **Empty TSV body** — the job died before producing output. Resubmit;
  if it persists, narrow the cell class.
- **5xx errors** — transient DDBJ-side issues. Retry with backoff.

## Primary citation

Zou Z. *et al.* (2024). ChIP-Atlas 3.0: a data-mining suite to explore
chromosome architecture together with large-scale regulome data.
*Nucleic Acids Research* 52(D1):D45–D53.
[doi:10.1093/nar/gkad884](https://doi.org/10.1093/nar/gkad884)
