---
name: okg-search
description: Query the statgen-analysis OKG for datasets matching a trait, EFO term,
  seed-dataset, or curated category, with controlled-vocab ancestry filtering (EUR/AFR/EAS/AMR/SAS/admixed/unknown),
  minimum-N power threshold, and sumstats-availability filter. Outputs a TSV ready
  to feed into `gwas-fetch GCST(X)` for batch chaining (e.g. genetic correlation with
  related immune traits). Auto-installs psycopg on first use and connects via the
  OKG_DSN env var; no MCP-server roundtrip overhead. Use when a user asks "what other
  GWAS exist for trait X", "what immune-trait studies are in the graph", "find sibling
  datasets like this one", or wants to bulk-pull a panel of GCSTs for rg / cross-trait
  analyses.
license: MIT
metadata:
  author: stephendorn
  version: '0.1'
  compatibility: Requires Python 3.9+, a reachable statgen-analysis Postgres (env
    var OKG_DSN; default postgres://postgres:okg@localhost:5449/statgen_analysis).
    Auto-installs `psycopg[binary]` on first run via uv pip / pip --user.
---

## DNAKit execution

Read [DNAKIT.md](DNAKIT.md) first for actual runtime, resource paths and validation status.


## What this skill does

Wraps four query modes against the statgen-analysis OKG and emits a TSV row per matching `dataset_metadata` node. Each row carries the GCST accession, study trait label, sample size, normalised ancestry codes, raw ancestry scope, genome build, release date, first author, PMID, and the GWAS Catalog sumstats URL (when present).

| Mode | Flag | What it does |
|---|---|---|
| Trait label | `--trait "inflammatory bowel"` | Fuzzy ILIKE match on `trait.attrs->>'label'`. Tries an exact match first, then substring. Returns datasets that `references` any matched trait. |
| EFO term | `--efo EFO:0003767` | Resolves a specific EFO term and returns every dataset whose trait `mapped_to` it. Catches BOTH the curated `trait:*` and auto-ingested `trait:gwas_catalog_scale:*` namespaces. |
| Seed dataset | `--seed-dataset dataset:gcst004131_ibd` | Follows `dataset --references--> trait --mapped_to--> EFO term`, then collects every dataset that shares the trait OR the EFO term. The "find sibling datasets" mode. |
| Curated category | `--category immune` | Filters the curated `trait:*` namespace by `attrs->>'category'` (e.g. `immune`, `psychiatric`, `metabolic`, `cardiovascular`). Skips the scale-ingest namespace because it doesn't have category populated. |

### Filters (applied after the mode resolution)

- `--ancestry EUR,EAS,...` — controlled-vocab codes `{EUR, AFR, EAS, AMR, SAS, admixed, unknown}`. The skill regex-parses each dataset's freeform `attrs->>'ancestry_scope'` into the set of codes that match (e.g. `"European, African"` → `{EUR, AFR}`). Multi-ancestry studies pass if ANY of their codes match the filter.
- `--min-n N` — keep only datasets with `n_cases + n_controls >= N`. Rows with no `n_cases`/`n_controls` populated are dropped by this filter (only the curated 5-trait subset has these fields today; scale-ingested datasets are mostly N=0).
- `--sumstats-only` — keep only datasets where `attrs->>'summary_statistics_url'` is non-empty. The other rows still have GCST accessions and you can usually fetch them via the GWAS Catalog FTP bucket, but the URL attribute means the OKG has verified the file location.

### Why direct Postgres, not MCP

Every other skill uses the MCP server (subprocess + stdio JSON-RPC). For bulk queries that hit hundreds of nodes (trait → dataset fan-outs), one MCP call per node = minutes of overhead. `okg-search` instead opens a Postgres connection via the `OKG_DSN` env var and runs a single JOIN — sub-second on the current ~205K-node graph. Same data, same access controls — Postgres is the source of truth that MCP wraps.

## EFO sibling expansion (default ON)

After resolving the initial trait set per mode, the skill walks **all seed traits → their EFO term(s) → every other trait that maps to those terms**. This bridges the two trait namespaces in the graph:

| Namespace | Carries `category` attr? | Populated by | Example |
|---|---|---|---|
| `trait:<curated_slug>` | yes (`immune`, `psychiatric`, ...) | hand-curated | `trait:crohns_disease`, `trait:asthma` |
| `trait:gwas_catalog_scale:<slug>` | no (default `gwas_catalog_mapped_trait`) | GWAS Catalog scale ingest (~thousands) | `trait:gwas_catalog_scale:crohn_s_disease` |

Without expansion, `--category immune` would only catch the curated namespace (5 traits → ~5 datasets, the curated set). With expansion, the same query bridges to scale-ingest siblings via the shared EFO terms: **5 curated immune traits → 9 expanded → 181 immune datasets**, of which 140 are EUR-ancestry. The same expansion improves `--trait` and `--efo` queries (catching multi-namespace duplicates with subtly different label phrasings).

Pass `--no-expand-efo` to keep only the seed traits (smaller, sharper result set; use when you specifically want the curated subset).

## Ancestry normalisation (regex; switches to graph attr when ready)

Today `dataset.attrs->>'ancestry_scope'` is freeform text. The skill ships a regex-based parser that emits the controlled vocabulary:

| Input substring (case-insensitive) | Code |
|---|---|
| european, caucasian, eur | EUR |
| african (not "african american"), afr | AFR |
| african american | AFR |
| east asian, chinese, japanese, korean, eas | EAS |
| south asian, indian, pakistani, bangladesh, sri lankan, sas | SAS |
| hispanic, latino, admixed american, amr | AMR |
| mixed, admixed (no qualifier), multi-ancestry, trans-ancestry | admixed |
| not reported, unknown, empty | unknown |

The OpenSpec change `normalize-dataset-ancestry-codes` (in `okg/openspec/changes/`) adds a derived `attrs->'ancestry_codes': List[str]` on every `dataset_metadata` node. Once that ships, the skill switches to a direct `attrs->'ancestry_codes' ? 'EUR'` JSONB containment check; the regex parser stays as the fallback for OKG generations that predate the change.

## Inputs you may need to elicit

One of:
- `--trait <label>` — exact-or-substring trait label match.
- `--efo <EFO:0000384>` — accepts `EFO:0000384`, `term:efo:0000384`, or bare `0000384`.
- `--seed-dataset <dataset:gcst...>` — OKG dataset node id (e.g. `dataset:gcst004131_ibd`).
- `--category <name>` — curated trait category (`immune`, `psychiatric`, `metabolic`, `cardiovascular`, etc.).

Optional filters:
- `--ancestry <codes>` — comma-separated controlled-vocab codes.
- `--min-n <int>` — minimum (n_cases + n_controls).
- `--sumstats-only` — STRICT: keep only datasets where the OKG has explicitly recorded a `summary_statistics_url` attribute. Today only the curated ~5 hand-imported datasets have this; the bulk 165K scale-ingest namespace doesn't. **Usually not the filter you want** — every GWAS Catalog accession is fetchable via `gwas-fetch` (which builds the URL from the GCST and the bucketed FTP path).
- `--has-gcst` — keep only rows with a non-empty GCST accession (effectively all `dataset_metadata` nodes). This is the **right filter when you plan to feed results into `gwas-fetch GCST<X>`**; it doesn't depend on the OKG having recorded an explicit `summary_statistics_url`.
- `--no-expand-efo` — disable the EFO sibling expansion (default ON; see "EFO expansion" below).
- `--limit <int>` — max rows in stdout TSV (default 50; the JSON manifest keeps all).
- `--json <path>` — write the full result set + query parameters to a JSON manifest.
- `--dsn <postgres://...>` — DSN override (default `$OKG_DSN` or `postgres://postgres:okg@localhost:5449/statgen_analysis`).

## How to execute

1. Run: `python3 scripts/okg_search.py <mode-flag> [filters]`.
2. The script writes a TSV to stdout (with the top `--limit` rows) and optionally a full JSON manifest to `--json <path>`.
3. Report back:
   - Trait nodes matched and dataset-count pre/post filtering (printed to stderr by the script).
   - Top results inline.
   - Suggested next step: pipe the GCST column into `gwas-fetch` for batch sumstats download (e.g. for genetic correlation via the `ldsc rg` subcommand).

## Refusal triggers

- No trait nodes matched the query → exit 1 with `[okg-search] no trait nodes matched`.
- OKG_DSN unreachable / wrong password → psycopg connection error propagates with the stack trace. Suggest checking the docker container (`docker ps | grep okg-pg-statgen-analysis`).
- Mode flags mutually exclusive: passing both `--trait` and `--efo` etc. is rejected by argparse.

## Examples

```bash
# 1. What other IBD datasets are in the graph?
python3 scripts/okg_search.py --trait "inflammatory bowel disease" --ancestry EUR

# 2. Find every dataset mapped to a specific EFO term
python3 scripts/okg_search.py --efo EFO:0003767 --sumstats-only --min-n 50000

# 3. "Find sibling datasets" — every study sharing IBD's trait/EFO with our seed
python3 scripts/okg_search.py \
    --seed-dataset dataset:gcst004131_ibd \
    --ancestry EUR --sumstats-only \
    --json results/ibd_siblings.json

# 4. All curated immune-category traits, EUR-ancestry only, with TSV out
python3 scripts/okg_search.py --category immune --ancestry EUR --limit 100

# 5. Bulk pipeline: find related EUR immune GWAS, pipe to gwas-fetch
python3 scripts/okg_search.py --seed-dataset dataset:gcst004131_ibd \
    --ancestry EUR --sumstats-only --min-n 20000 --limit 25 \
  | tail -n +2 | cut -f2 \
  | while read gcst; do
        python3 ~/.claude/skills/gwas-fetch/scripts/gwas_fetch.py --gcst "$gcst"
    done
```

## Output schema

TSV columns (in order): `node_id`, `gcst`, `trait`, `n_total`, `ancestry_codes`, `ancestry_scope`, `build`, `release_date`, `author`, `pmid`, `sumstats_url`.

`ancestry_codes` is a comma-separated list (e.g. `EUR,AFR`); `ancestry_scope` is the raw freeform string the codes were parsed from.

## JSON manifest schema

`<--json path>`:

```json
{
  "skill": "okg-search",
  "version": "0.1",
  "query": {
    "trait": null, "efo": null,
    "seed_dataset": "dataset:gcst004131_ibd", "category": null,
    "ancestry": ["EUR"], "min_n": 20000,
    "sumstats_only": true, "limit": 25
  },
  "n_results": 5,
  "results": [
    {
      "node_id": "dataset:gcst004131_ibd", "gcst": "GCST004131",
      "trait": "Inflammatory bowel disease",
      "ancestry_scope": "european",
      "n_total": 59957, "n_cases": "25042", "n_controls": "34915",
      "sumstats_url": "https://www.ebi.ac.uk/gwas/rest/api/studies/GCST004131/associations",
      "release_date": "2017-01-09T00:00:00Z",
      "build": "GRCh37",
      "author": "...", "pmid": "28067908",
      "matched_traits": ["trait:gwas_catalog_scale:inflammatory_bowel_disease"],
      "ancestry_codes": ["EUR"]
    }
  ],
  "captured_at": "2026-05-21T..."
}
```

## Notes and edge cases

- **Two trait namespaces.** The graph has both `trait:<curated_slug>` (manually curated, ~few hundred, populated `category`) and `trait:gwas_catalog_scale:<slug>` (auto-ingested, ~thousands, no `category`). Modes `--seed-dataset` and `--efo` walk the EFO term to merge both; modes `--trait` and `--category` may surface one but not the other.
- **EFO subtree traversal is not yet supported.** The OKG ingests EFO terms as flat nodes without `is_a` / `subclass_of` parent edges, so `--efo-subtree EFO:0005140` (autoimmune disease → all child traits) can't be answered today. Filed as a future OpenSpec change.
- **Scale-ingested datasets often have N=0.** The gen-37 (2026-04-27) GWAS Catalog snapshot populates `accession` + `release_date` + `disease_trait` + `ancestry_scope` but not `n_cases`/`n_controls` for the bulk `dataset:gwas_catalog_scale:*` namespace. The curated 5 datasets (IBD, CD, UC, asthma, Crohn's) DO have them. `--min-n` therefore filters out most scale-ingest hits; remove it if you want a wider net.
- **Provenance.** The skill does not emit per-row OKG provenance (the rows ARE the OKG node IDs). When the result feeds downstream skills (`gwas-fetch`, `ldsc rg`, etc.), those skills record their own provenance citing the same node IDs.
