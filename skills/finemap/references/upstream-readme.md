# statgen-skills

A collection of Claude Code / agent skills for statistical-genetics workflows, built to the [Agent Skills](https://agentskills.io) standard. Each skill is a self-contained directory with a `SKILL.md` + helper scripts.

## Available skills

| Skill | Purpose |
|---|---|
| [`gwas-fetch/`](https://github.com/svdorn/statgen-skills/blob/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/gwas-fetch) | Fetch harmonised GWAS Catalog sumstats by accession (`GCST*`). OKG-first metadata resolution with REST-API fallback. Caches downloads + writes a `.fetch.json` provenance manifest. Compose with `liftover` for build alignment. |
| [`liftover/`](https://github.com/svdorn/statgen-skills/blob/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/liftover) | KG-aware genomic liftover for GWAS sumstats. Resolves source/target build via the OKG when available; falls back to the GWAS Catalog REST API. Refuses unknown-provenance lifts and emits a coverage-gap proposal stub. Supports `--gcst <accession>` for direct GWAS-Catalog fetch+lift. |
| [`okg-search/`](https://github.com/svdorn/statgen-skills/blob/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/okg-search) | Query the statgen-analysis OKG for datasets matching a trait, EFO term, seed-dataset, or curated category. Filters: controlled-vocab ancestry (`EUR/AFR/EAS/AMR/SAS/admixed`), min-N power threshold, sumstats availability. EFO sibling expansion bridges the curated `trait:*` and scale-ingest `trait:gwas_catalog_scale:*` namespaces. TSV output is ready to feed into `gwas-fetch GCST<X>` for batch chaining (e.g. `rg` panels). Connects via direct Postgres (`OKG_DSN`); auto-installs `psycopg[binary]` on first use. |
| [`manhattan-qq/`](https://github.com/svdorn/statgen-skills/blob/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/manhattan-qq) | Quick QC for a GWAS sumstats file: Manhattan + QQ plots + a summary JSON with genomic-inflation λ_GC, GWS + suggestive hit counts, and a list of independent lead SNPs (greedy ±window-kb pruning). Auto-detects common column-name conventions (harmonised `hm_*`, LDSC munged, bare `CHR/BP/P/SNP`). Subsamples low-significance variants for fast plotting on full HM3- or 1KG-scale sumstats. Writes a `.manhattan_qq.json` sidecar citing the input dataset's OKG node. |
| [`ldsc/`](https://github.com/svdorn/statgen-skills/blob/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/ldsc) | Run LD Score Regression: `munge` raw sumstats, estimate SNP heritability (`h2`), compute genetic correlation (`rg`), and stratified / cell-type-specific heritability (`partitioned`, wrapping LDSC `--h2-cts` for the Finucane 2018 multi-tissue / GenoSkyline-Plus / any `.ldcts` panel). Auto-installs the CBIIT LDSC fork; auto-detects `--signed-sumstats` from the file header (hm_beta → 0, hm_odds_ratio → 1, etc.); auto-supplies `--merge-alleles` from the on-disk HM3 snplist so indels don't crash downstream `rg`; writes a `.ldsc.json` sidecar with OKG provenance and parsed key results. |
| [`finemap/`](https://github.com/svdorn/statgen-skills/blob/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/finemap) | Fine-mapping with multiple interchangeable backends: **sushie** (`mancusolab/sushie`, single- and multi-ancestry SuSiE / SuShiE) as the default; **SuSiE 2.0** (`StatFunGen/susieR`, Wang 2025) as the default SuSiE-inf implementation; and **FinucaneLab/fine-mapping-inf** (Cui 2024) as an opt-in (`--susie-inf-impl finucane`). Five subcommands: `susie` / `sushie` (individual genotypes), `region` (GWAS sumstats + reference VCF/PLINK; auto-detects rsid-vs-chrpos SNP-ID format of the panel), `sumstats` (precomputed Z + LD with method dispatch), and `recommend` (OKG-driven method choice given a dataset_metadata node). Each backend auto-installs on first use; writes a `.finemap.json` sidecar citing the method / software / paper / ld_panel / dataset OKG nodes for whichever backend ran. |
| [`prs/`](https://github.com/svdorn/statgen-skills/blob/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/prs) | Polygenic risk scoring via a pluggable `--method` flag. Currently implements **SBayesRC** (Zheng 2024) through the `zhilizheng/SBayesRC` R package. OKG-resolves the matching HapMap3 LD eigendecomposition by `--ancestry eur\|eas\|afr` (or explicit `--okg-ld-panel-id`), auto-pulls `N` from `--okg-dataset-id`, converts harmonised / LDSC / GWAS-SSF sumstats to COJO format, **backfills missing effect-allele frequencies from the LD reference's `snp.info`** (GWAS-Catalog harmonised files often ship all-NA EAF; without backfill SBayesRC's `tidy()` drops every row), runs SBayesRC's 3-step `tidy → impute → sbayesrc` pipeline, and writes per-SNP PRS weights + `.prs.json` sidecar citing 5 OKG nodes (method, software, paper, ld_panel, dataset). LDpred2 / PRS-CS easily wired in behind the same flag. |
| [`twas/`](https://github.com/svdorn/statgen-skills/blob/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/twas) | Transcriptome-wide association study via **TWAS-FUSION** (Gusev et al. 2016, *Nat Genet*) against the GTEx v8 multi-tissue weights. Resolves the (`--tissue`, `--ancestry`) pair to an OKG `dataset:fusion_gtex_v8_<eur\|all>:<tissue>` node, downloads + caches the panel + the FUSION 1000G EUR LDREF on first use, runs `FUSION.assoc_test.R` per chromosome, concatenates per-gene results, and writes a `.twas.json` sidecar citing the method / software / paper / ld_panel / tissue / cohort / dataset OKG nodes. Auto-installs FUSION's R deps (`plink2R`, `glmnet`, `here`, `coloc`) on first use. With `--okg-trait-id`, queries the OKG's `tissue → trait relevant_to` mapping and iterates over every tissue relevant to that trait. |
| [`locuszoom/`](https://github.com/svdorn/statgen-skills/blob/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/locuszoom) | LocusZoom-style regional association plots via the [`locuszoomr`](https://github.com/myles-lewis/locuszoomr) R package. Selects a window by gene symbol, region (`chr:start-end`), or lead-SNP rsid; fetches pairwise r² from the **LDlink** REST API (1000G, configurable population); renders a multi-track PDF + PNG with the association scatter on top and Ensembl gene tracks below. With `--finemap-pip` or `--finemap-sidecar`, adds a middle panel of per-variant PIPs coloured by credible-set membership — the canonical 3-panel fine-mapped-locus figure. Requires a free LDlink API token — the skill prompts on first use and caches it under `~/.cache/locuszoom/ldlink_token`. Writes a `.locuszoom.json` sidecar with locus selection, LD coverage stats, and OKG provenance (propagated from the upstream `.finemap.json` when supplied). |
| [`variant-annotate/`](variant-annotate/) | Per-variant functional + tissue annotation against gnomAD (allele frequencies, predicted consequence, gene context) and GTEx (top eQTL / sQTL hits). Used by `finemap region` to enrich the printed credible-set summary table with rsid → consequence / gene / AF / eQTL annotations. Caches REST responses under `~/.cache/variant-annotate/`. |

## Installation

The skills follow the [Agent Skills](https://agentskills.io/specification) specification, so they should work with any compatible agent (Claude Code, Cursor, Goose, OpenCode, etc.).

### Claude Code (user-level)

Symlink each skill into `~/.claude/skills/`:

```bash
git clone https://github.com/svdorn/statgen-skills.git
cd statgen-skills
./install.sh                  # symlinks all skills into ~/.claude/skills/
```

Or symlink a single skill:

```bash
ln -s ~/statgen-skills/liftover ~/.claude/skills/liftover
```

After symlinking, restart Claude Code (or the agent client) so it picks up the new skills at startup.

### Other agentskills.io-compatible clients

Point your client at `statgen-skills/<skill>` or copy individual skills into the client's skills directory. Each skill's `SKILL.md` documents its `compatibility` field for any environment requirements.

## Skill development

Each skill is structured per the agentskills.io spec:

```
<skill-name>/
├── SKILL.md          # required: YAML frontmatter (name, description, ...) + body
├── scripts/          # executable helpers (Python, bash, etc.)
├── references/       # detailed docs loaded on demand
└── assets/           # templates, schemas, etc.
```

Validate a skill before committing:

```bash
# optional: uses agentskills' reference validator
skills-ref validate ./liftover
```

## OKG integration

Several skills here are designed around an **OKG** (Operational Knowledge Graph) — the upstream reference is [github.com/mitdbg/okg](https://github.com/mitdbg/okg/tree/main/deployments/statgen-analysis), but any fork or compatible deployment works. To use OKG-aware features, set the `OKG_REPO` environment variable to your local clone:

```bash
git clone https://github.com/mitdbg/okg ~/Lab/KG/okg   # or your fork
export OKG_REPO=~/Lab/KG/okg
```

The skills then spawn the OKG MCP server at `$OKG_REPO/deployments/<deployment>/server.py` (default deployment: `statgen-analysis`) and query for provenance metadata before acting. When the OKG can't resolve a needed entity, skills emit a coverage-gap proposal stub instead of guessing — see each skill's `references/COVERAGE_GAPS.md` for the workflow.

If `$OKG_REPO` is not set and a skill is invoked with an `--okg-*` flag, the skill refuses with a clear error pointing you here. Skills can still run **without** the OKG by passing all required metadata explicitly (e.g. `--source hg38 --target hg19` for liftover).

## License

MIT — see [LICENSE](https://github.com/svdorn/statgen-skills/blob/dc5f0c60d2995aebfe58b3f8a07b46f27562f4e0/LICENSE).
