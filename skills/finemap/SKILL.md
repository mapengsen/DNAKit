---
name: finemap
description: Fine-map GWAS loci with SuShiE/SuSiE or SuSiE-inf using individual genotypes
  or summary statistics and matched LD panels. Supports single/multiple ancestry,
  credible sets, optional gnomAD/GTEx annotation and provenance manifests. Use for
  posterior inclusion probabilities, credible sets or fine-mapping method selection;
  solver/reference installation is separate.
license: MIT
metadata:
  author: stephendorn
  version: '0.3'
  source: agentskills.io/specification
  sushie_repo: https://github.com/mancusolab/sushie
  sushie_docs: https://mancusolab.github.io/sushie/
  sushie_paper_doi: 10.1038/s41588-025-02262-7
  susie_inf_default_impl: susie2
  susie_inf_susie2_repo: https://github.com/StatFunGen/susieR
  susie_inf_susie2_label: SuSiE 2.0 (Wang 2025)
  susie_inf_finucane_repo: https://github.com/FinucaneLab/fine-mapping-inf
  susie_inf_finucane_paper_doi: 10.1038/s41588-023-01597-3
  compatibility: Requires Python 3.9+, git, and ~500 MB of disk for the sushie clone
    + tutorial data. Auto-installs sushie via `pip install .` (or `uv pip install`).
    `--method susie_inf` defaults to installing R (brew/apt) + StatFunGen/susieR (SuSiE
    2.0). `--susie-inf-impl finucane` opts into FinucaneLab/fine-mapping-inf (Python).
    OKG provenance optional, gated on $OKG_REPO.
---

## DNAKit execution

Read [DNAKIT.md](DNAKIT.md) first for actual runtime, resource paths and validation status.


## What this skill does

Wraps the `sushie finemap` CLI behind a thin orchestrator that:

1. **Installs sushie** on first use (clones `mancusolab/sushie` to `~/.cache/sushie/repo`, runs `pip install .`).
2. **Verifies the install** by running the bundled 3-ancestry tutorial against `<repo>/data/`. Once-only; results recorded under `~/.cache/sushie/verify_install/`.
3. **Runs fine-mapping** in one of four modes:
   - **`susie`**: single-population, individual-level (one VCF + one phenotype file).
   - **`sushie`**: multi-population, individual-level (K VCFs + K phenotype files).
   - **`region`**: single-locus from GWAS sumstats + reference genotypes. Pass `--gwas-sumstats`, `--chrom/--start/--end`, `--N`, and one of `--ref-vcf` / `--ref-plink` / `--ref-bgen`; sushie computes the LD matrix internally. This is the simplest mode when you have a GWAS sumstats file and an ancestry-matched genotype reference panel (e.g. 1000G EUR).
   - **`sumstats`**: precomputed Z-score files + precomputed LD matrices per ancestry (K=1 SuSiE-RSS, K≥2 SuShiE) via `infer_sushie_ss` — use when you've already computed LD yourself.
4. **Records provenance** in a `.finemap.json` sidecar with sushie version (git SHA), bundled CLI flags, parsed CS / PIP summary, and OKG node IDs (method, software, paper, ld_panel, dataset where applicable).

### Choosing between `region` and `sumstats`

| You have… | Use `region` | Use `sumstats` |
|---|---|---|
| GWAS sumstats + reference VCF / PLINK / BGEN | ✓ (sushie computes LD) | requires precomputed LD |
| Z-scores + precomputed LD `.npy` per locus | n/a | ✓ |
| Multiple ancestries simultaneously (SuShiE) | not yet (use `sumstats`) | ✓ (pass K files per flag) |

For LD-source options + how to get a 1000G EUR (hg19/hg38) reference VCF, see [references/LD_OPTIONS.md](references/LD_OPTIONS.md).

### SNP-ID matching (the #1 "no common SNPs" cause)

Sushie joins the GWAS sumstats to the reference panel by **exact string equality on the SNP ID** — so the GWAS-side ID must use the same convention as the reference panel's `.bim` / VCF ID column. Reference panels in the wild use one of two conventions:

| Panel | SNP-ID format | Example |
|---|---|---|
| HM3 / 1000G PLINK files / on-disk `~/.cache/ldsc/.../1000G.EUR.QC.<chr>.bim` | **rsid** | `rs55768441` |
| 1000G Phase3 VCFs from FTP / sushie's `--ref-1000g` slice | **`.` (no ID)** → forced chr:pos | `1:66742081` |
| UKB / custom panels | varies | either |

The `region` subcommand has a `--snp-id-format {auto, rsid, chrpos}` flag (default `auto`):

- **`auto`**: peeks at the first 500 variant IDs in the reference panel's `.bim` / VCF. If >50% match `rs<digits>`, emits rsid-style GWAS IDs; otherwise emits `chr:pos`.
- **`rsid`** / **`chrpos`**: forces one format regardless of what the panel looks like (useful when the auto-detection is wrong, e.g. a panel where IDs are a mix of rsids and chr:pos).

Failure mode: if you see `ValueError: Ancestry 1: No common SNPs between GWAS and genotype data`, it's almost always an ID-format mismatch — try `--snp-id-format` override. The skill prints the chosen format on every run for diagnosis.

### Reference-panel sources for `region`

The `--ref-plink` / `--ref-vcf` / `--ref-bgen` flags accept any ancestry-matched genotype panel. Reuseable on-disk sources you may already have:

- **`~/.cache/ldsc/ldsc_inputs/for_enrichment/genotype/1000G.EUR.QC.<chr>.{bed,bim,fam}`** — 1000G Phase 3 EUR PLINK, all chroms, rsid-keyed. Pass as `--ref-plink .../1000G.EUR.QC.<chr>` (one chromosome at a time for `region`). This is the path of least resistance for EUR fine-mapping when the ldsc skill cache is populated.
- **`--ref-1000g {eur,eas,afr,sas,amr}`** — convenience: pulls the 1000G Phase 3 VCF slice for your `--chrom:--start-end` from the public FTP, subsets to the requested superpopulation, and uses it as the LD reference. First-use installs htslib + bcftools and seeds `~/.cache/finemap/1000g_<build>_<pop>/`. The slice ships with `ID=.` so the skill auto-detects `chrpos` format.
- **Your own custom PLINK/VCF**: any panel works as long as it's on the same genome build as the GWAS and has compatible SNP-ID conventions.

## First-run setup (verify-install)

The first time the skill is invoked (or any time the user passes `--verify-install`), the skill runs:

```bash
cd ~/.cache/sushie/repo/data/
sushie finemap \
    --pheno EUR.pheno AFR.pheno EAS.pheno \
    --vcf vcf/EUR.vcf vcf/AFR.vcf vcf/EAS.vcf \
    --covar EUR.covar AFR.covar EAS.covar \
    --output ~/.cache/sushie/verify_install/test_result
```

This is the canonical tutorial from https://mancusolab.github.io/sushie/manual.html — 3 ancestries (489 EUR, 639 AFR, 481 EAS), 123/129/113 HapMap SNPs.

Success criteria:
- `test_result.cs.tsv` exists and contains at least one credible set
- `test_result.weight.tsv` exists
- No tracebacks in stderr

A `~/.cache/sushie/verify_install/.verified` sentinel file is written on success; future invocations skip the verification unless `--re-verify` is passed.

## Inputs you may need to elicit

For `susie` (single-ancestry):
- `--vcf <path>` — genotype VCF
- `--pheno <path>` — phenotype file (tab-separated, sushie format)
- `--covar <path>` (optional) — covariates
- `--out <prefix>` — output prefix

For `sushie` (multi-ancestry, K ≥ 2):
- `--vcf <eur.vcf> <afr.vcf> [...]` — K VCFs (order matters, must match pheno order)
- `--pheno <eur.pheno> <afr.pheno> [...]` — K phenotype files
- `--covar <eur.covar> <afr.covar> [...]` (optional)
- `--out <prefix>`

For `sumstats` (GWAS sumstats + LD; method-dispatched):
- `--z <z1.tsv> [<z2.tsv> ...]` — z-score TSV per ancestry (or BETA+SE columns; computed on the fly)
- `--ld <ld1.npy> [<ld2.npy> ...]` — LD matrix per ancestry (`.npy` or whitespace TSV, square, same SNP order as `--z`)
- `--n <n1> [<n2> ...]` — sample size per ancestry
- `--out <prefix>`
- `--method {sushie, susie_inf, susie_ash, auto}` — default `sushie`. `susie_inf` and `susie_ash` are K=1 only.
- `--okg-dataset-id <node>` (optional) — required for `--method auto`; recorded in the manifest regardless of method.
- See [references/SUMMARY_STATS.md](references/SUMMARY_STATS.md) for file formats and examples.

### Choosing a method (`sumstats` mode)

| `--method` | When | Multi-ancestry? | Backend impl |
|---|---|---|---|
| `sushie` (default) | K=1 standard GWAS, OR K≥2 multi-ancestry | ✓ K≥2 native | mancusolab/sushie (Python, py310 venv) |
| `susie_inf` | Moderate-large N single-ancestry GWAS (N≥100k); want to model an infinitesimal polygenic background | K=1 only | Picked by `--susie-inf-impl` (default `susie2` = StatFunGen/susieR; opt-in `finucane` = FinucaneLab/fine-mapping-inf) |
| `auto` | "let the OKG pick" — needs `--okg-dataset-id` | depends on graph signals | n/a |

`--susie-inf-impl` (only relevant when `--method susie_inf`):

| Value | Backend | Paper / OKG citation |
|---|---|---|
| `susie2` (**default**) | StatFunGen/susieR — SuSiE 2.0 `susie_rss(..., unmappable_effects="inf")` | `paper:susie_ash_2025` (Wang 2025) + `software:susie_ash` |
| `finucane` | FinucaneLab/fine-mapping-inf — Cui 2024 Python | `paper:susie_inf_2024` + `software:susie_inf` |

Both impls implement the same SuSiE-inf statistical model (`method:susie_inf` in the OKG). SuSiE 2.0 converges more cleanly on real GWAS in our IL23R benchmark (recovers the canonical R381Q variant rs11209026 with 2 CSs; FinucaneLab+MLE missed it). FinucaneLab is preserved for regression-testing / continuity with prior runs.

Method dispatch by `--method auto` (or the standalone `recommend` subcommand):

```
K=2+ inputs              -> sushie  (multi-ancestry SuShiE)
OKG dataset has multi-ancestry siblings -> sushie  (suggest user supply per-ancestry sumstats)
n_total >= 100_000       -> susie_inf  (default impl: susie2)
otherwise                -> sushie  (K=1, SuSiE-RSS-equivalent)
```

For `recommend` (OKG-driven method choice; no actual fine-mapping):
- `--okg-dataset-id <node>` — the GWAS dataset (e.g. `dataset:gcst90245992_height`)
- `--n-inputs N` (optional, default 1) — how many ancestry inputs you plan to supply
- Prints `{method, rationale, signals}` JSON.

## How to execute

1. **Run the script**: `python3 scripts/finemap.py <subcommand> <flags>`.
2. **First-run install + verify** are automatic on the first invocation; subsequent runs skip them.
3. **Report back**:
   - Output prefix and files generated (`<out>.cs.tsv`, `<out>.weight.tsv`, `<out>.cv.tsv`, `<out>.log`)
   - Number of credible sets, mean CS size, top-PIP variant per CS
   - For `sushie`: the estimated cross-ancestry effect-size correlation matrix from `<out>.corr.tsv`
   - For `region`: a printed **rsID-enriched credible-set table** (CS index, PIP, rsID, chr:pos, effect/non-effect alleles, β, SE, p) — joined from sushie's chr:pos-keyed CS output back to the original GWAS sumstats via the `<out>.gwas_extras.tsv` sidecar the orchestrator writes during the COJO-style conversion. One row per CS member, sorted by CS then descending PIP within CS. Markdown format so it's directly pasteable into reports / chat.
   - Path to `.finemap.json` sidecar with full provenance

## Examples

```bash
# Verify install (run the 3-ancestry tutorial once)
python3 scripts/finemap.py --verify-install

# Single-ancestry SuSiE on one locus
python3 scripts/finemap.py susie \
    --vcf GWAS/genotypes/MSH3.vcf \
    --pheno GWAS/phenos/T01.pheno \
    --out results/finemap/MSH3_T01

# Multi-ancestry SuShiE across EUR + AFR + EAS
python3 scripts/finemap.py sushie \
    --vcf GWAS/genotypes/EUR_MSH3.vcf GWAS/genotypes/AFR_MSH3.vcf GWAS/genotypes/EAS_MSH3.vcf \
    --pheno GWAS/phenos/EUR_T01.pheno GWAS/phenos/AFR_T01.pheno GWAS/phenos/EAS_T01.pheno \
    --out results/finemap/MSH3_T01_3anc

# Sumstats fine-mapping (single ancestry, SuSiE-RSS)
python3 scripts/finemap.py sumstats \
    --z  GWAS/sumstats/MSH3.z.tsv \
    --ld LD/MSH3.ld.npy \
    --n  361194 \
    --out results/finemap/MSH3_sumstats

# Sumstats fine-mapping (multi-ancestry SuShiE)
python3 scripts/finemap.py sumstats \
    --z  EUR_MSH3.z.tsv  AFR_MSH3.z.tsv  EAS_MSH3.z.tsv \
    --ld EUR_MSH3.ld.npy AFR_MSH3.ld.npy EAS_MSH3.ld.npy \
    --n  361194           48000           92000 \
    --out results/finemap/MSH3_3anc_sumstats

# Sumstats fine-mapping with SuSiE-inf (default impl: SuSiE 2.0 / StatFunGen, K=1)
python3 scripts/finemap.py sumstats \
    --method susie_inf \
    --z  GWAS/sumstats/MSH3.z.tsv \
    --ld LD/MSH3.ld.tsv \
    --n  1597374 \
    --out results/finemap/MSH3_susie_inf

# Same model via the original FinucaneLab impl (regression-testing / continuity)
python3 scripts/finemap.py sumstats \
    --method susie_inf --susie-inf-impl finucane \
    --z  GWAS/sumstats/MSH3.z.tsv \
    --ld LD/MSH3.ld.npy \
    --n  139274 \
    --out results/finemap/MSH3_susie_inf_finucane

# Auto: let the OKG pick the method given a dataset_metadata node
OKG_REPO=~/Lab/KG/okg python3 scripts/finemap.py sumstats \
    --method auto \
    --okg-dataset-id dataset:gcst90245992_height \
    --z  GWAS/sumstats/height.z.tsv \
    --ld LD/height.ld.npy \
    --n  1597374 \
    --out results/finemap/height_auto

# Standalone method recommendation (no run, just the choice + rationale)
OKG_REPO=~/Lab/KG/okg python3 scripts/finemap.py recommend \
    --okg-dataset-id dataset:gcst004988_brca

# Pin to a specific sushie commit (reproducibility)
python3 scripts/finemap.py sushie --sushie-commit <sha> --vcf ... --pheno ... --out ...
```

## Sidecar manifest schema

`<out>.finemap.json`:

```json
{
  "subcommand": "sushie",
  "output_prefix": "results/finemap/MSH3_T01_3anc",
  "sushie_repo": "https://github.com/mancusolab/sushie",
  "sushie_commit": "<git sha>",
  "n_ancestries": 3,
  "inputs": {
    "vcf":   ["EUR.vcf", "AFR.vcf", "EAS.vcf"],
    "pheno": ["EUR.pheno", "AFR.pheno", "EAS.pheno"],
    "covar": ["EUR.covar", "AFR.covar", "EAS.covar"]
  },
  "summary": {
    "n_cs": 3,
    "mean_cs_size": 5.7,
    "max_pip": 0.98,
    "cross_ancestry_correlation": [[1.0, 0.85, 0.74], [0.85, 1.0, 0.79], [0.74, 0.79, 1.0]]
  },
  "okg_node_ids": {
    "method":   "method:sushie",
    "software": "software:sushie",
    "paper":    "paper:sushie_2025"
  },
  "captured_at": "2026-05-13T..."
}
```

For `susie` mode, `okg_node_ids.method` is `method:susie_finemapping` (the canonical SuSiE method), with `software:sushie` still as the operational implementation — the same dual-provenance pattern used by the `ldsc` skill.

## Notes and edge cases

- **sushie's CLI is `sushie finemap`**; this skill wraps it with `python3 finemap.py {susie,sushie}` for the individual-level path, and calls `sushie.infer_ss.infer_sushie_ss` directly for the `sumstats` path.
- **Summary-statistics inputs**: pass per-ancestry z-score TSVs, LD matrices, and sample sizes via `--z`, `--ld`, `--n`. See [references/SUMMARY_STATS.md](references/SUMMARY_STATS.md) for file formats, an alignment guide, and notes on LD-panel pairing.
- **VCF format**: sushie accepts plink/VCF/BGEN. The skill's flags use `--vcf` for simplicity; pass `--plink <prefix>` to use plink format instead.
- **Single-ancestry mode** is sushie's degenerate K=1 case; results are equivalent to SuSiE-RSS at numerical tolerance.

For details on sushie's bundled test data and the verification step, see [references/EXAMPLES.md](references/EXAMPLES.md).
