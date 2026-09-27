---
name: ldsc
description: Run LD Score Regression (LDSC) for sumstats munging, SNP heritability
  (h2), genetic correlation (rg), and stratified/cell-type-specific partitioned heritability
  (partitioned, Finucane 2015/2018 S-LDSC via `--h2-cts`). Auto-installs LDSC from
  the CBIIT Python 3 / Mac-compatible fork and caches the canonical EUR LD-score reference
  files on first use. Resolves method/software/paper provenance from the OKG when
  $OKG_REPO is set; otherwise runs without OKG provenance. Use when a user asks to
  munge GWAS sumstats, estimate heritability from sumstats, compute genetic correlation
  between two GWAS traits, or rank tissues / cell types by heritability enrichment
  (S-LDSC / Finucane 2018 multi-tissue).
license: MIT
metadata:
  author: stephendorn
  version: '0.1'
  source: agentskills.io/specification
  ldsc_fork: https://github.com/CBIIT/ldsc
  ldsc_upstream: https://github.com/bulik/ldsc
  compatibility: Requires Python 3.9+, git, and an internet connection on first run
    to clone CBIIT/ldsc and download the EUR LD-score reference (~250 MB extracted).
    Auto-installs the LDSC Python deps (numpy, pandas, scipy, bitarray) into the active
    environment via uv pip / pip --user.
---

## DNAKit execution

Read [DNAKIT.md](DNAKIT.md) first for actual runtime, resource paths and validation status.


## What this skill does

Wraps the three core LDSC workflows behind a single CLI:

| Subcommand | Wraps | Purpose |
|---|---|---|
| `munge` | `munge_sumstats.py` | Standardize a raw GWAS TSV into LDSC's `.sumstats.gz` format. |
| `h2` | `ldsc.py --h2` | Estimate SNP heritability from one set of munged sumstats. |
| `rg` | `ldsc.py --rg` | Estimate genetic correlation between two sets of munged sumstats. |
| `partitioned` | `ldsc.py --h2-cts` | Stratified / cell-type-specific heritability (S-LDSC). Joint-fits a baseline (Finucane 2015 baseline or Finucane 2018 baseline-LD v2.2) with a per-tissue annotation panel listed in a `.ldcts` file (Finucane 2018 multi-tissue gene expression / chromatin, GenoSkyline-Plus, etc.). One LDSC invocation produces a per-tissue ranked `.cell_type_results.txt`. |

On first use the skill:
1. Clones `https://github.com/CBIIT/ldsc` to `~/.cache/ldsc/repo/` (override via `--repo-cache`).
2. Installs LDSC's Python deps into the active env (`numpy`, `pandas`, `scipy`, `bitarray`).
3. Downloads `eur_w_ld_chr.tar.bz2` from the Broad and extracts to `~/.cache/ldsc/ld_scores/eur_w_ld_chr/` (override via `--ld-scores-dir`).

Subsequent runs reuse the cache; no network needed unless the user passes `--refresh`.

## Inputs you may need to elicit

Ask via `AskUserQuestion` when missing:

- **`munge`**: input TSV path, output prefix, and (if not standard) the sumstats column names (`--snp`, `--a1`, `--a2`, `--p`, `--n`, etc. mirror LDSC's `munge_sumstats.py` flags). **The signed-sumstats column is auto-detected** from the file header — `hm_beta` → `hm_beta,0`, `hm_odds_ratio` → `hm_odds_ratio,1`, else falls back to bare `beta` / `odds_ratio`. Pass `--signed-sumstats <col>,<null>` to override. **`--merge-alleles` is also auto-supplied** when an HM3 snplist is found at one of `~/.cache/ldsc/ldsc_inputs/w_hm3.snplist`, `~/.cache/ldsc/w_hm3.snplist`, or `~/.cache/ldsc/ld_scores/w_hm3.snplist`. This restricts the output to biallelic HM3 SNPs (LDSC's canonical convention for h² / rg) and prevents the silent `KeyError: '<indel-allele>'` crash that LDSC's `rg` throws when indels with multi-character allele strings survive munge. Override by passing `--merge-alleles <path>` after `--` in the extras.
- **`h2`**: input `.sumstats.gz` and (optionally) a partitioned LD-score reference if the user wants stratified h2. Default reference is `eur_w_ld_chr/`.
- **`rg`**: two `.sumstats.gz` files (`--in1`, `--in2`) and the LD-score reference.
- **`partitioned`**: input `.sumstats.gz`, plus three required reference paths:
  - `--baseline-prefix` — per-chrom baseline ldscore prefix, e.g. `~/.cache/ldsc/ldsc_inputs/finucane_2018/baselineLD_v2.2/baselineLD.` (trailing dot intended; LDSC appends `<chr>.l2.ldscore.gz`).
  - `--ldcts` — path to a `.ldcts` file (shipped inside `Multi_tissue_*.tgz`, e.g. `Multi_tissue_gene_expr.ldcts`) listing tissues and their per-chrom ldscore prefixes.
  - `--weights-prefix` — per-chrom HM3 regression-weight prefix, e.g. `~/.cache/ldsc/ldsc_inputs/for_enrichment/weights/weights.hm3_noMHC.` (the on-disk `weights.hm3_noMHC` is equivalent to the Finucane 2018 `weights_hm3_no_hla` — same MHC/HLA region excluded).
  - Optional `--samp-prev` / `--pop-prev` for liability-scale h² (same semantics as `h2`).

## How to execute

1. **Run the script**: `python3 scripts/ldsc.py <subcommand> <flags>`. The script handles install, cache, deps, sidecar manifest.

2. **Report back**:
   - Output paths (`.log`, `.sumstats.gz`, `.results`, etc. per subcommand)
   - Key numbers parsed from the LDSC `.log`:
     - munge: number of variants retained / dropped, mean chi-squared
     - h2: SNP h² estimate ± SE, intercept ± SE, ratio
     - rg: rg estimate ± SE, p-value, h² of each trait
   - Path to `<output>.ldsc.json` sidecar manifest with OKG provenance + reference-file SHA-256s + LDSC version (git commit).

3. **OKG provenance**:
   - If `$OKG_REPO` set: query `get_node("method:ldsc")` and `get_node("software:ldsc")` at the pinned generation; embed the node IDs in the sidecar manifest.
   - If `$OKG_REPO` unset: the skill still runs; manifest just omits `okg_node_ids`.

4. **Refusal triggers**:
   - munge: input has fewer than 100,000 variants after standardization (LDSC needs ≥200k for stable estimates; flag a warning at 200k, refuse below 100k).
   - h2: mean chi-squared < 1.02 (the GWAS is unpowered for h² estimation — flag, but don't refuse).
   - rg: either input's h² < 0 or SE > estimate (the input is too noisy for genetic correlation).

## Examples

```bash
# 1. Munge raw sumstats into LDSC format
python3 scripts/ldsc.py munge \
    --in GWAS/raw/T01.tsv.gz \
    --out GWAS/munged/T01 \
    --N 419013   # or auto-read from a column

# 2. Estimate heritability
python3 scripts/ldsc.py h2 \
    --in GWAS/munged/T01.sumstats.gz \
    --out h2/T01

# 3. Genetic correlation between two traits
python3 scripts/ldsc.py rg \
    --in1 GWAS/munged/T01.sumstats.gz \
    --in2 GWAS/munged/T02.sumstats.gz \
    --out rg/T01_T02

# 4. With OKG provenance (records method:ldsc, software:ldsc, paper:ldsc_2015 in the manifest)
OKG_REPO=~/Lab/KG/okg python3 scripts/ldsc.py h2 --in T01.sumstats.gz --out h2/T01

# 5. Cell-type-specific S-LDSC (Finucane 2018 multi-tissue gene expression panel)
python3 scripts/ldsc.py partitioned \
    --in ldsc/ibd.sumstats.gz \
    --out sldsc/ibd_multitissue \
    --baseline-prefix ~/.cache/ldsc/ldsc_inputs/finucane_2018/baselineLD_v2.2/baselineLD. \
    --ldcts        ~/.cache/ldsc/ldsc_inputs/finucane_2018/Multi_tissue_gene_expr_1000Gv3_ldscores/Multi_tissue_gene_expr.ldcts \
    --weights-prefix ~/.cache/ldsc/ldsc_inputs/for_enrichment/weights/weights.hm3_noMHC. \
    -- --samp-prev 0.4177 --pop-prev 0.005
```

**Outputs from `partitioned`:**
- `<out>.cell_type_results.txt` — TSV with `Name Coefficient Coefficient_std_error Coefficient_P_value` per tissue (LDSC's native output).
- `<out>.log` — LDSC log including overall h² + intercept fit.
- `<out>.ldsc.json` — sidecar manifest with the top tissues, the `.ldcts` and baseline/weights prefixes used, OKG provenance, and (when available) liability-scale settings.

## Stratified LDSC (S-LDSC) — what it does and how to run it

### What S-LDSC actually does

Plain LDSC (`h2`) fits **one** regression slope of χ²(SNP) on the SNP's LD-score (the sum of r² with all SNPs in the reference panel). The slope is proportional to per-SNP heritability, and the intercept catches inflation that is *not* polygenic (stratification, cryptic relatedness, model misspecification).

S-LDSC (`partitioned`, via LDSC's `--h2-cts`) generalises that regression to **many overlapping annotations**: each SNP gets a vector of per-annotation LD-scores, and the regression returns a per-annotation coefficient τ_c. τ_c > 0 means SNPs in annotation c carry more per-SNP h² than the baseline expectation. Two common questions it answers:

1. **Functional partition.** Joint-fit ~97 baseline-LD v2.2 categories (coding, enhancer, conserved, etc.) against a trait's sumstats; ask which functional category carries the most h² per SNP. Output: `.results` table with `Enrichment` (= Prop_h2 / Prop_SNPs) and `Coefficient` per category.
2. **Tissue / cell-type partition (Finucane 2018 multi-tissue, this skill's `partitioned` subcommand).** For each tissue in a `.ldcts` file, joint-fit `[baseline categories] + [tissue-specific annotation] + [control annotation]` and test whether the tissue annotation's coefficient is positive. Output: `.cell_type_results.txt` with one row per tissue, columns `Name Coefficient Coefficient_std_error Coefficient_P_value`.

The `partitioned` subcommand drives (2). For (1), pass `--ref-ld-chr` and `--overlap-annot` through `--extra` to `h2`; the per-category enrichment ends up in `<out>.results`.

### When to run S-LDSC

- **Pick the right TWAS / coloc tissue from the data, not from prior intuition.** S-LDSC ranks tissues by h² enrichment in one regression; the top GTEx tissue feeds directly into TWAS-FUSION via the `twas` skill.
- **Sanity-check trait biology.** A trait that doesn't enrich anywhere sensible (e.g. height enriching only in brain) signals a sumstats or model issue.
- **Compare two traits' tissue profiles** by running the same `--ldcts` against both munged sumstats and lining up the rankings.

Do *not* run S-LDSC when mean χ² < 1.1 — coefficients become too noisy to rank. Use plain `h2` first to check power.

### Required reference data (Finucane 2018 multi-tissue panel)

The canonical files used to live at `alkesgroup.broadinstitute.org/LDSCORE/`, which is **deprovisioned as of 2024-2025**. Current locations:

| File | Purpose | Source |
|---|---|---|
| `1000G_Phase3_baselineLD_v2.2_ldscores.tgz` (~676 MB) | 97-category joint-fit functional baseline | Zenodo DOI `10.5281/zenodo.10515792` OR `gs://broad-alkesgroup-public-requester-pays/LDSCORE/` |
| `Multi_tissue_gene_expr_1000Gv3_ldscores.tgz` (~2.2 GB) | 53 GTEx v6p + 152 Franke-lab tissues, plus the canonical `Multi_tissue_gene_expr.ldcts` | `gs://broad-alkesgroup-public-requester-pays/LDSCORE/LDSC_SEG_ldscores/` (requester-pays, no public mirror) |
| `weights_hm3_no_hla` (~50 MB) | HM3 regression weights with the MHC/HLA region excluded | same requester-pays bucket — OR reuse on-disk `weights.hm3_noMHC.<chr>` (equivalent for this purpose) |

Requester-pays = needs a GCP project + `gcloud auth login` + `gsutil -u <project-id> cp gs://...`. Cost ~$0.18 for 1.5 GB (egress at $0.12/GB). Setup ~10 min from scratch.

Recommended on-disk layout:

```
~/.cache/ldsc/ldsc_inputs/finucane_2018/
├── baselineLD_v2.2/                              # 22 chroms × {annot.gz, l2.ldscore.gz, l2.M, l2.M_5_50}
├── Multi_tissue_gene_expr_1000Gv3_ldscores/      # 205 tissues × 22 chroms
└── Multi_tissue_gene_expr.ldcts                  # tissue → ldscore-prefix mapping
```

### Running `partitioned`

```bash
python3 scripts/ldsc.py partitioned \
    --in ldsc/<trait>.sumstats.gz \
    --out sldsc/<trait>_multitissue \
    --baseline-prefix ~/.cache/ldsc/ldsc_inputs/finucane_2018/baselineLD_v2.2/baselineLD. \
    --ldcts          ~/.cache/ldsc/ldsc_inputs/finucane_2018/Multi_tissue_gene_expr.ldcts \
    --weights-prefix ~/.cache/ldsc/ldsc_inputs/for_enrichment/weights/weights.hm3_noMHC. \
    -- --samp-prev <case_frac> --pop-prev <population_prevalence>
```

Notes on the args:
- `--baseline-prefix` / `--weights-prefix` end in a trailing `.` — LDSC appends `<chr>.l2.ldscore.gz` to that prefix.
- `--ldcts` is the per-tissue mapping file. The skill rewrites it internally with absolute paths so LDSC can resolve the per-tissue prefixes regardless of cwd (upstream LDSC has a long-standing bug here).
- The trailing `--` separates skill flags from arbitrary `ldsc.py` flags (passed through verbatim). Use it for liability-scale conversion (`--samp-prev <case-fraction> --pop-prev <population-prevalence>`).
- Runtime: one LDSC invocation but it loops over tissues internally. The full 205-tissue Finucane 2018 panel takes **~15–20 min** on a modern laptop for a typical HM3-munged sumstats. Each tissue runs an independent S-LDSC fit; they don't share state, so this can't be trivially parallelised inside LDSC.

### Interpreting the output

`<out>.cell_type_results.txt` is a TSV with one row per tissue (LDSC's native format):

| Column | Meaning |
|---|---|
| `Name` | Tissue label from the `.ldcts` (matches GTEx-FUSION tissue names for the 53 GTEx tissues) |
| `Coefficient` | τ_c — per-SNP h² contribution of the tissue annotation above what the baseline-LD model + control annotation explain |
| `Coefficient_std_error` | Block-jackknife SE for the coefficient |
| `Coefficient_P_value` | One-sided p (H₀: τ_c ≤ 0). Small p = tissue annotation is *enriched* for h². |

**Bonferroni correction.** With 205 tissues in the Finucane 2018 panel, the Bonferroni-significance threshold is `0.05 / 205 = 2.44e-04`. Apply per `.ldcts` (don't combine multiple `.ldcts` panels into one multiple-testing budget without thought).

**Top hits feed TWAS.** The top GTEx tissue (filter `Name` to those without the Franke `A##.###` prefix) maps 1:1 to the `--tissue` flag of the `twas` skill — e.g., `Whole_Blood`, `Spleen`, `Colon_Sigmoid`. Pass that to TWAS-FUSION as the expression-weight panel.

**Caveats:**
- The coefficient is *not* an enrichment ratio. A small p doesn't mean the tissue carries most of the h² — it means the per-SNP contribution from this tissue annotation is significantly positive *after* controlling for the baseline. Use the per-category `.results` from a baseline-LD-only run to get enrichment ratios.
- Highly correlated tissue annotations (e.g. GTEx Brain regions, or GTEx `Whole_Blood` vs Franke `Blood`) will share the signal — both can come up significant for an immune trait because they have ~50% overlap, not because they're independently informative. Treat ranking as a coarse guide, not a fine-grained ordering.
- S-LDSC needs ≥200K SNPs in the regression after merging with HM3 weights. Sumstats with smaller overlap will throw a warning; LDSC still runs, but coefficients get noisy.

## Sidecar manifest schema

Every run writes `<out_prefix>.ldsc.json` with:

```json
{
  "subcommand": "h2",
  "input": "...", "output_prefix": "...",
  "ldsc_repo": "https://github.com/CBIIT/ldsc",
  "ldsc_commit": "<git sha>",
  "ld_scores_dir": "~/.cache/ldsc/ld_scores/eur_w_ld_chr",
  "ld_scores_sha256": "<sha256 of the tarball>",
  "key_results": {
    "h2": 0.123, "h2_se": 0.014,
    "intercept": 1.02, "intercept_se": 0.008,
    "ratio": 0.05, "mean_chi2": 1.34
  },
  "okg_node_ids": {"method": "method:ldsc", "software": "software:ldsc", "paper": "paper:ldsc_2015"},
  "captured_at": "2026-05-13T..."
}
```

## Notes and edge cases

- **CBIIT fork vs bulik/ldsc**: the OKG currently records `bulik/ldsc` as the canonical software. The CBIIT fork is functionally equivalent for the workflows here but Python-3 / macOS-arm64 compatible. Override via `--repo-url <url>` if you want the upstream or another fork.
- **LD-score files** are downloaded once from the Broad's public mirror. The default (EUR) covers EUR-ancestry GWAS; for trans-ancestry / non-EUR use, pass `--ld-scores-dir` to a panel matching your sumstats' ancestry. The skill flags ancestry-mismatch warnings when the input's `population` column doesn't match the LD-scores' assumed ancestry.
- **Allele alignment** during munge: LDSC handles allele flips automatically when an effect-allele frequency column is present. The skill prints how many SNPs were ambiguous and dropped.
- **Reproducibility**: pin LDSC by passing `--repo-commit <sha>` to clone-checkout a specific revision. Default = `main` HEAD.

See [references/LDSC_REFS.md](references/LDSC_REFS.md) for the canonical EUR/EAS/AFR LD-score-file URLs and citation guidance.
