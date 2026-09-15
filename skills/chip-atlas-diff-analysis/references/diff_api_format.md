# ChIP-Atlas Diff Analysis — API Specification

This document records the request/response contract used by the WABI Diff
Analysis endpoint. The skill's scripts talk to this endpoint directly; this
file is the reference for when something breaks or you need to adapt the
call from another script.

## Endpoint

```
POST https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/
```

The same host serves the Peak Enrichment Analysis endpoint. The
`antigenClass` field is what selects "diff" mode: pass `diffbind` or `dmr`
instead of an antigen class name.

## POST parameters

Every submission must include all of these as form fields. Empty strings
are required for the placeholder fields — do not omit them.

| Field | Required value or example | Notes |
|-------|--------------------------|-------|
| `address` | `""` | Reserved field; must be present but empty |
| `format` | `"text"` | **Diff analysis returns text** (enrichment uses JSON) |
| `result` | `"www"` | Reserved field |
| `genome` | `hg38` / `hg19` / `mm10` / `mm9` / `rn6` / `dm6` / `dm3` / `ce11` / `ce10` / `sacCer3` | One of ten supported assemblies |
| `antigenClass` | `"diffbind"` or `"dmr"` | Selects pipeline (edgeR vs metilene) |
| `cellClass` | `"empty"` | Fixed; no cell-type filter applies in diff mode |
| `threshold` | `"50"` | Fixed for diff mode |
| `typeA` | `"srx"` | Tells the API that `bedAFile` contains experiment IDs |
| `bedAFile` | Newline-separated IDs | At least 2 of `SRX…`, `ERX…`, `DRX…`, or `GSM…` |
| `typeB` | `"srx"` | Same as `typeA` |
| `bedBFile` | Newline-separated IDs | At least 2 of `SRX…`, `ERX…`, `DRX…`, or `GSM…` |
| `permTime` | `"1"` | Fixed |
| `title` | Free text | Job title shown in the log file |
| `descriptionA` | Free text | Label for group A in plots / IGV |
| `descriptionB` | Free text | Label for group B in plots / IGV |

## Supported assemblies

| Organism | Assemblies |
|----------|-----------|
| *Homo sapiens* | hg38, hg19 |
| *Mus musculus* | mm10, mm9 |
| *Rattus norvegicus* | rn6 |
| *Drosophila melanogaster* | dm6, dm3 |
| *Caenorhabditis elegans* | ce11, ce10 |
| *Saccharomyces cerevisiae* | sacCer3 |

## Accepted experiment ID formats

The API expects DDBJ/EBI/NCBI run accessions (any of `SRX`, `ERX`, `DRX`)
or GEO sample accessions (`GSM`). Each list is collapsed into a
newline-separated string when sent over HTTP.

Minimum two IDs per group; three or more is strongly recommended so edgeR
can estimate dispersion. Mixing accession types within one group is
allowed.

## Response format — initial submission

The body is tab/colon-separated lines. Example:

```
requestId	wabi_chipatlas_2025-0517-1942-12-987-654321
parameters	null
current-time	2025-05-17 19:42:12
start-time
current-state
```

Pull `requestId` out of the first line by splitting on whitespace or
colon. The skill's `query_chip_atlas_diff.submit_diff_job` accepts either
`requestId: …` or `requestId\t…` formatting.

## Polling for completion

```
GET https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/{requestId}?info=status
```

The body looks like:

```
status: running
```

or

```
status: finished
```

Failed jobs return `status: error` or `status: failed`. Diff jobs commonly
take 2–10 minutes; allow up to 15 for large groups. Poll at 10–30 second
intervals.

## Retrieving the results

```
GET https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/{requestId}?info=result&format=zip
```

The response body is a ZIP archive. Files live inside a `wabi_result/`
subdirectory:

| Member | Format | Purpose |
|--------|--------|---------|
| `wabi_result.bed` | 8-col TSV (no header) | Main result — see `output_format.md` |
| `wabi_result.igv.bed` | BED9 + URL-encoded GFF3 | IGV visualization; orange/blue color encodes direction |
| `wabi_result.xml` | XML (IGV session) | Loads BigWig tracks + the BED in IGV |
| `wabi_result.log` | text | Sample list, pipeline parameters, runtime |
| `README.md` | markdown | Short description of the files |

## Differences vs the Enrichment endpoint

The diff endpoint and the enrichment endpoint share a host but differ in
ways that matter for the client:

1. Diff uses `format=text` on submit; enrichment uses `format=json`.
2. Diff fixes `cellClass=empty` — there is no cell-type filter.
3. `bedAFile` / `bedBFile` carry **experiment IDs**, not gene symbols or
   regions.
4. Results come back as a **ZIP** archive of files, not a TSV. Retrieve
   the binary content and unzip into a working directory.
5. Diff jobs run slower than enrichment (minutes, not seconds).

## Failure modes worth handling

- `HTTP 400` — almost always a malformed parameter. Reread the table
  above; common mistakes are missing `address` / `result` (the API
  rejects calls without them) or a `genome` outside the supported set.
- `HTTP 5xx` — transient DDBJ outage. Back off and retry.
- The polling endpoint sometimes returns an HTML error page when the
  job is still being scheduled; the script ignores intermittent failures
  for up to the configured timeout.
- A successful poll followed by an invalid ZIP usually means the job
  finished but its result file is still being staged. Wait a few seconds
  and retry the GET.

## Primary citation

Zou Z. *et al.* (2024). ChIP-Atlas 3.0: a data-mining suite to explore
chromosome architecture together with large-scale regulome data.
*Nucleic Acids Research* 52(D1):D45–D53.
[doi:10.1093/nar/gkad884](https://doi.org/10.1093/nar/gkad884)
