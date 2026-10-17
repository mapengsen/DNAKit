# GWAS

This category contains **29** documented capabilities. All 29 compute results natively in DNAKit and none returns `not_implemented` in a normal call: summary-statistics handling (`read`, `standardize`, `validate`, `filter`, `harmonize`, `write`, `export`), QC and plot data (`calculate`, `apply_genomic_control`, `compute_qq_data`, `compute_manhattan_data`, `compute_locus_plot_data`, `plot`), locus definition and fine mapping (`select_lead_snps`, `merge_nearby_significant_hits`, `define_association_loci`, `annotate_lead_snps`, `build_credible_set`, `conditional_gwas`, `joint_gwas_analysis`), meta-analysis (`fixed_effect_meta_analysis`, `random_effect_meta_analysis`), heritability and genetic correlation (`estimate_snp_heritability`, `estimate_genetic_correlation`), polygenic scores (`clump_and_threshold_prs`, `evaluate_prs`, `calibrate_prs_risk`), and the native association scan and simulation (`run`, `simulate_gwas_data`). numpy is optional (it only speeds up the genomic relationship matrix). The one exception is `run` with `backend` set to an external program (PLINK, REGENIE, GEMMA, GCTA): it returns `status: "not_implemented"` because DNAKit does not drive external programs. The functions are exported through the DNAKit namespace; functions with a compatible JSON call plan are also exposed through the Agent/MCP registry.

## 1) annotate_lead_snps

- **Purpose:** Annotate genome-wide-significant lead SNPs with the nearest and overlapping genes from a user-supplied gene table.

- **Method:** Lead SNPs are either given (`leads`) or selected as in `select_lead_snps` (`p_threshold`, `window`). For each lead the distance to every gene on the same chromosome is 0 inside a gene and otherwise the distance to the nearest gene boundary; the overlapping genes, the nearest gene (only if within `max_distance`) and the genes within `flank` bp are reported. No database is queried and no gene model is bundled.

- **API:** `dnakit.gwas.annotate_lead_snps(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. `options`: `genes` (required; records with `name`, `chrom`, `start`, `end` and optional `strand`), `leads` (explicit variant ids; default selected with `p_threshold` = 5e-8 and `window` = 500000), `flank` (100000), `max_distance` (1000000).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `annotations` (per lead: the standardised record fields, `neglog10p`, `n_nearby_significant`, `overlapping_genes`, `nearest_gene`, `distance_to_nearest_gene`, `relation`, `genes_within_flank`), `n_leads` and `method`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import annotate_lead_snps

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
genes = [
    {"name": "GENE1", "chrom": "1", "start": 90000, "end": 120000, "strand": "+"},
    {"name": "GENE2", "chrom": "1", "start": 880000, "end": 940000, "strand": "-"},
    {"name": "GENE3", "chrom": "2", "start": 500000, "end": 520000, "strand": "+"},
]
result = annotate_lead_snps(rows, options={"genes": genes, "flank": 100000})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "annotate_lead_snps",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "annotations": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.12,
        "se": 0.02,
        "p": 2e-09,
        "eaf": 0.3,
        "neglog10p": 8.698970004336019,
        "n_nearby_significant": 0,
        "overlapping_genes": [
          "GENE1"
        ],
        "nearest_gene": "GENE1",
        "distance_to_nearest_gene": 0,
        "relation": "within_gene",
        "genes_within_flank": [
          "GENE1"
        ]
      },
      {
        "variant_id": "rs3",
        "chrom": "1",
        "pos": 900000,
        "effect_allele": "G",
        "other_allele": "A",
        "beta": -0.15,
        "se": 0.025,
        "p": 2e-09,
        "eaf": 0.4,
        "neglog10p": 8.698970004336019,
        "n_nearby_significant": 0,
        "overlapping_genes": [
          "GENE2"
        ],
        "nearest_gene": "GENE2",
        "distance_to_nearest_gene": 0,
        "relation": "within_gene",
        "genes_within_flank": [
          "GENE2"
        ]
      }
    ],
    "n_leads": 2,
    "method": "nearest gene by interval distance from the supplied gene table"
  }
}
```

## 2) apply_genomic_control

- **Purpose:** Apply genomic control: divide chi-square statistics by the inflation factor lambda_GC and recompute p-values.

- **Method:** lambda_GC = median(chi2)/0.4549364, where chi2 is the 1-df chi-square statistic implied by each p-value and 0.4549364 is the median of the 1-df chi-square distribution. The used factor is lambda_GC, `options['lambda_gc']` if given, and with `only_inflated=True` (default) max(1, lambda) so that nothing is deflated. Each record gets `p_gc` (from chi2/lambda), `se_gc` = se*sqrt(lambda) and `z_gc` = z/sqrt(lambda).

- **API:** `dnakit.gwas.apply_genomic_control(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. It needs `p` (or `z`, or `beta` with `se`). `options`: `lambda_gc` (overrides the estimate), `only_inflated` (True).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `lambda_gc` (the estimate), `lambda_used`, `records` (the standardised records with `p_gc`, `se_gc`, `z_gc` added) and `n_records`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import apply_genomic_control

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = apply_genomic_control(rows, options={"lambda_gc": 1.2})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "gwas",
  "operation": "apply_genomic_control",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "lambda_gc": 26.65064698571636,
    "lambda_used": 1.2,
    "records": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.12,
        "se": 0.02,
        "p": 2e-09,
        "eaf": 0.3,
        "n": 5000,
        "z": 6.0,
        "p_gc": 4.3695935095779474e-08,
        "se_gc": 0.02190890230020664,
        "z_gc": 5.477225575051661
      },
      {
        "variant_id": "rs2",
        "chrom": "1",
        "pos": 150000,
        "effect_allele": "C",
        "other_allele": "T",
        "beta": 0.09,
        "se": 0.02,
        "p": 6.8e-06,
        "eaf": 0.25,
        "n": 5000,
        "z": 4.5,
        "p_gc": 3.994693946558859e-05,
        "se_gc": 0.02190890230020664,
        "z_gc": 4.107919181288747
      },
      "..."
    ],
    "n_records": 6
  }
}
```

## 3) build_credible_set

- **Purpose:** Build a fine-mapping credible set from per-variant beta and se using Wakefield approximate Bayes factors (single causal variant).

- **Method:** log ABF = 0.5*(ln(V/(V+W)) + W/(V+W)*z^2) with V = se^2, z = beta/se and W = prior_sd^2. With a uniform prior the posterior inclusion probability is PIP = ABF/sum(ABF); the credible set is the smallest set of top-PIP variants whose PIPs sum to at least `coverage`. It assumes exactly one causal variant in the region and ignores LD.

- **API:** `dnakit.gwas.build_credible_set(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. The table must hold the variants of ONE locus and include `beta` and `se`. `options`: `prior_sd` (0.15), `coverage` (0.95), `region` (`chrom:start-end`, to subset the table).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `credible_set` (variant ids), `set_size`, `cumulative_pip`, `coverage_target`, `prior_sd`, `variants` (every variant with `log_abf` and `pip`, sorted by PIP) and `method`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import build_credible_set

locus = [
    {"variant_id": "v1", "chrom": "1", "pos": 1000, "beta": 0.30, "se": 0.05},
    {"variant_id": "v2", "chrom": "1", "pos": 1100, "beta": 0.28, "se": 0.05},
    {"variant_id": "v3", "chrom": "1", "pos": 1200, "beta": 0.10, "se": 0.05},
    {"variant_id": "v4", "chrom": "1", "pos": 1300, "beta": 0.02, "se": 0.05},
    {"variant_id": "v5", "chrom": "1", "pos": 1400, "beta": -0.05, "se": 0.05},
]
result = build_credible_set(locus, options={"prior_sd": 0.15, "coverage": 0.95})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "gwas",
  "operation": "build_credible_set",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 5,
    "text_length": 374,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9edc672e99be6afc5d0c531ed01d24ab9a6f80d910fe1821d4838015578a5c21"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "credible_set": [
      "v1",
      "v2"
    ],
    "set_size": 2,
    "cumulative_pip": 0.9999992874121282,
    "coverage_target": 0.95,
    "prior_sd": 0.15,
    "variants": [
      {
        "variant_id": "v1",
        "chrom": "1",
        "pos": 1000,
        "log_abf": 15.048707453502972,
        "pip": 0.8897307255825493
      },
      {
        "variant_id": "v2",
        "chrom": "1",
        "pos": 1100,
        "log_abf": 12.960707453502978,
        "pip": 0.11026856182957885
      },
      "..."
    ],
    "method": "Wakefield approximate Bayes factors, single causal variant per region, uniform prior; log ABF = 0.5*(ln(V/(V+W)) + W/(V+W)*z^2) with V = se^2, W = prior_sd^2"
  }
}
```

## 4) calculate

- **Purpose:** Calculate a summary statistic of GWAS summary statistics; the default is a table summary including the genomic inflation factor lambda_GC.

- **Method:** `summary`: number of variants and of p-values, minimum p and its variant, lambda_GC = median(chi2)/0.4549364, mean and median chi2 (chi2 implied by p), number of variants with p <= `significance` and number of chromosomes. `lambda_gc` and `n_significant` give just those values; `z_scores`, `neglog10p` (-log10 p) and `maf` (min(eaf, 1-eaf)) give one value per variant; `adjusted_p` gives the Bonferroni-adjusted p (p x m, capped at 1) and Benjamini-Hochberg q-values.

- **API:** `dnakit.gwas.calculate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. `options`: `statistic` = `summary` (default), `lambda_gc`, `z_scores`, `neglog10p`, `maf`, `n_significant` or `adjusted_p`; `significance` (5e-8).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `statistic` and: for `summary`, `n_variants`, `n_with_p`, `min_p`, `lead_variant`, `lambda_gc`, `mean_chi2`, `median_chi2`, `significance_threshold`, `n_significant`, `n_chromosomes`; for `lambda_gc`, `value` and `median_chi2`; for `n_significant`, `threshold` and `value`; for the per-variant statistics, `values` (and `n_tests` for `adjusted_p`). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import calculate

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = calculate(rows, options={"statistic": "summary"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "calculate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "statistic": "summary",
    "n_variants": 6,
    "n_with_p": 6,
    "min_p": 2e-09,
    "lead_variant": "rs1",
    "lambda_gc": 26.65064698571636,
    "mean_chi2": 16.24053827366549,
    "median_chi2": 12.124350013504214,
    "significance_threshold": 5e-08,
    "n_significant": 2,
    "n_chromosomes": 2
  }
}
```

## 5) calibrate_prs_risk

- **Purpose:** Calibrate a polygenic score to absolute disease risk with a logistic recalibration model.

- **Method:** Fits risk = expit(a + b*z) by logistic regression of the 0/1 phenotype on the standardised score z (individuals with a missing value are dropped; at least 10 are needed). With `prevalence` the intercept is shifted by logit(prevalence) - logit(sample prevalence) to undo case-control sampling. Reports the Brier score and a calibration table of mean predicted versus observed risk in `n_bins` score-ranked bins.

- **API:** `dnakit.gwas.calibrate_prs_risk(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'score': [...], 'phenotype': [0/1, ...]}`. `options`: `prevalence` (population prevalence in (0,1)), `n_bins` (10), `new_scores` (list of scores to convert to risk).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `intercept`, `slope_per_sd`, `score_mean`, `score_sd`, `sample_prevalence`, `prevalence_adjusted`, `predicted_risk`, `calibration_table` (`bin`, `n`, `mean_predicted`, `observed`), `brier_score`, `mean_predicted`, `new_score_risk` (with `new_scores`) and `method`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import calibrate_prs_risk

score = [-1.2, -0.8, -0.5, -0.3, 0.0, 0.1, 0.3, 0.6, 0.9, 1.4, -1.0, -0.2, 0.2, 0.7, 1.1, 1.6]
phenotype = [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1]
result = calibrate_prs_risk(
    {"score": score, "phenotype": phenotype},
    options={"n_bins": 4, "new_scores": [-1.0, 0.0, 1.0]},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "gwas",
  "operation": "calibrate_prs_risk",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 160,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "bf1d3cbdd5f0dd2548bfac81709c02e63edd6840e4f64c0e340776a98b2507e6"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "intercept": -0.46097257162881095,
    "slope_per_sd": 2.360136952448344,
    "score_mean": 0.18125,
    "score_sd": 0.835239486614468,
    "sample_prevalence": 0.4375,
    "prevalence_adjusted": false,
    "predicted_risk": [
      0.012567774392076724,
      0.037917350382972455,
      0.08424795828956246,
      "..."
    ],
    "calibration_table": [
      {
        "bin": 1,
        "n": 4,
        "mean_predicted": 0.039159842928464444,
        "observed": 0.0
      },
      {
        "bin": 2,
        "n": 4,
        "mean_predicted": 0.23107124190436826,
        "observed": 0.5
      },
      "..."
    ],
    "brier_score": 0.13582393605542953,
    "mean_predicted": 0.4375,
    "new_score_risk": [
      0.021906288649246153,
      0.27425689620543464,
      0.8644285709560183
    ],
    "method": "logistic recalibration risk = expit(a + b * standardised score); with `prevalence`, a is shifted by logit(prevalence) - logit(sample prevalence)"
  }
}
```

## 6) clump_and_threshold_prs

- **Purpose:** Compute clumping-and-thresholding (C+T) polygenic scores from summary statistics and dosage genotypes.

- **Method:** PLINK-style greedy clumping: variants are taken in ascending p; each becomes an index variant and absorbs the not-yet-clumped variants on the same chromosome within `clump_window` bp whose r^2 with it, computed from the LD reference genotypes, is at least `clump_r2`. For each p-value threshold the score of a sample is sum(beta*dosage) over the index variants with p <= threshold, a missing dosage being replaced by the variant mean. The LD reference defaults to the target genotypes themselves.

- **API:** `dnakit.gwas.clump_and_threshold_prs(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'summary_statistics': table with variant_id, chrom, pos, beta, p, 'genotypes': samples x variants effect-allele dosages in [0, 2], 'sample_ids': [...], 'variant_ids': [...]}`. `options`: `p_thresholds` (`[5e-8, 1e-5, 1e-3, 0.05, 0.5, 1]`), `clump_r2` (0.1), `clump_window` (250000), `ld_reference` (a genotype mapping; default the target genotypes).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `clumps` (`index`, `p`, `beta`, `members`), `thresholds` (per threshold: `p_threshold`, `n_variants`, `variants` and `scores` {sample_id: score}), `clump_r2`, `clump_window` and `method`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import clump_and_threshold_prs

summary = [
    {"variant_id": "v1", "chrom": "1", "pos": 1000, "beta": 0.30, "p": 1e-9},
    {"variant_id": "v2", "chrom": "1", "pos": 2000, "beta": 0.25, "p": 1e-7},
    {"variant_id": "v3", "chrom": "1", "pos": 3000, "beta": -0.20, "p": 1e-4},
]
genotypes = [[0, 0, 2], [1, 1, 1], [2, 2, 0], [1, 1, 0], [0, 0, 1], [2, 2, 2]]
result = clump_and_threshold_prs(
    {"summary_statistics": summary, "genotypes": genotypes,
     "sample_ids": ["a", "b", "c", "d", "e", "f"], "variant_ids": ["v1", "v2", "v3"]},
    options={"p_thresholds": [1e-8, 1e-3], "clump_r2": 0.5},
)
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "gwas",
  "operation": "clump_and_threshold_prs",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 411,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3c2cc5467f66341cafe87e4b0959583a1cf91dd85b0704c13b015e61777908d7"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "clumps": [
      {
        "index": "v1",
        "p": 1e-09,
        "beta": 0.3,
        "members": [
          "v2"
        ]
      },
      {
        "index": "v3",
        "p": 0.0001,
        "beta": -0.2,
        "members": []
      }
    ],
    "thresholds": [
      {
        "p_threshold": 1e-08,
        "n_variants": 1,
        "variants": [
          "v1"
        ],
        "scores": {
          "a": 0.0,
          "b": 0.3,
          "c": 0.6,
          "d": 0.3,
          "e": 0.0,
          "f": 0.6
        }
      },
      {
        "p_threshold": 0.001,
        "n_variants": 2,
        "variants": [
          "v1",
          "v3"
        ],
        "scores": {
          "a": -0.4,
          "b": 0.09999999999999998,
          "c": 0.6,
          "d": 0.3,
          "e": -0.2,
          "f": 0.19999999999999996
        }
      }
    ],
    "clump_r2": 0.5,
    "clump_window": 250000,
    "method": "PLINK-style greedy clumping by p-value (r^2 from the LD reference, same chromosome, within window bp), then score = sum(..."
  }
}
```

## 7) compute_locus_plot_data

- **Purpose:** Compute the points of a regional (LocusZoom-style) association plot around a lead variant.

- **Method:** The region is the lead position +/- `flank` bp (start clamped at 1) on the lead's chromosome; every variant in it becomes a point with -log10 p. If an LD matrix or reference genotypes are supplied, r^2 with the lead is computed (squared Pearson correlation of dosages, or of the given r) and each point is put in a band (0-0.2, 0.2-0.4, 0.4-0.6, 0.6-0.8, 0.8-1, `lead`, or `unknown` without LD information). No LD is fabricated.

- **API:** `dnakit.gwas.compute_locus_plot_data(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. `options`: `lead` (variant id; default the smallest p), `flank` (250000), `ld` (`{'variant_ids': [...], 'matrix': [[...]]}` correlations r) or `genotypes` (reference dosages) to colour by r^2.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `lead` (`variant_id`, `chrom`, `pos`, `p`), `region` (`chrom`, `start`, `end`) and `points` (`variant_id`, `pos`, `neglog10p`, `r2`, `ld_band`). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import compute_locus_plot_data

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = compute_locus_plot_data(rows, options={"lead": "rs1", "flank": 100000})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "compute_locus_plot_data",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "lead": {
      "variant_id": "rs1",
      "chrom": "1",
      "pos": 100000,
      "p": 2e-09
    },
    "region": {
      "chrom": "1",
      "start": 1,
      "end": 200000
    },
    "points": [
      {
        "variant_id": "rs1",
        "pos": 100000,
        "neglog10p": 8.698970004336019,
        "r2": 1.0,
        "ld_band": "lead"
      },
      {
        "variant_id": "rs2",
        "pos": 150000,
        "neglog10p": 5.167491087293763,
        "r2": null,
        "ld_band": "unknown"
      }
    ]
  }
}
```

## 8) compute_manhattan_data

- **Purpose:** Compute Manhattan-plot coordinates (-log10 p and a cumulative genomic x position) from summary statistics.

- **Method:** x = chromosome offset + pos, where the offset of a chromosome is the sum of the maximum observed positions of the preceding chromosomes (chromosomes ordered 1..22, X, Y, MT). y = -log10 p. Points weaker than `min_neglog10p` are dropped; above `max_points` the most significant half of the allowed points is kept exactly and the rest is thinned. The genome-wide (5e-8) and suggestive (1e-5) thresholds and the chromosome label centres are reported.

- **API:** `dnakit.gwas.compute_manhattan_data(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. It needs `chrom`, `pos` and `p`. `options`: `min_neglog10p` (0), `max_points` (200000).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `points` (`variant_id`, `chrom`, `pos`, `x`, `neglog10p`), `chromosomes` (`chrom`, `offset`, `length`, `center`), `total_length`, `truncated`, `thresholds` (`genome_wide`, `suggestive`) and `x_definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import compute_manhattan_data

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = compute_manhattan_data(rows)
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "gwas",
  "operation": "compute_manhattan_data",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "points": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100000,
        "x": 100000,
        "neglog10p": 8.698970004336019
      },
      {
        "variant_id": "rs2",
        "chrom": "1",
        "pos": 150000,
        "x": 150000,
        "neglog10p": 5.167491087293763
      },
      "..."
    ],
    "chromosomes": [
      {
        "chrom": "1",
        "offset": 0,
        "length": 900000,
        "center": 450000.0
      },
      {
        "chrom": "2",
        "offset": 900000,
        "length": 90000,
        "center": 945000.0
      }
    ],
    "total_length": 990000,
    "truncated": false,
    "thresholds": {
      "genome_wide": {
        "p": 5e-08,
        "neglog10p": 7.301029995663981
      },
      "suggestive": {
        "p": 1e-05,
        "neglog10p": 5.0
      }
    },
    "x_definition": "cumulative base-pair position: chromosome offset (sum of preceding chromosomes' maximum observed positions) + pos"
  }
}
```

## 9) compute_qq_data

- **Purpose:** Compute quantile-quantile plot data (observed versus expected -log10 p) and lambda_GC.

- **Method:** For p sorted ascending with rank 1 the smallest, expected = -log10((rank - 0.5)/n) and observed = -log10 p. Above `max_points` the most significant half of the points is kept exactly and the rest is thinned. lambda_GC = median(chi2)/0.4549364.

- **API:** `dnakit.gwas.compute_qq_data(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. It needs `p`. `options`: `max_points` (50000).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `n`, `points` (`rank`, `expected`, `observed`), `lambda_gc`, `downsampled` and `expected_definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import compute_qq_data

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = compute_qq_data(rows)
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "gwas",
  "operation": "compute_qq_data",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n": 6,
    "points": [
      {
        "rank": 6,
        "expected": 0.0377885608893998,
        "observed": 0.2076083105017461
      },
      {
        "rank": 5,
        "expected": 0.12493873660829993,
        "observed": 0.49894073778224846
      },
      "..."
    ],
    "lambda_gc": 26.65064698571636,
    "downsampled": false,
    "expected_definition": "-log10((rank - 0.5) / n) with rank 1 = smallest p"
  }
}
```

## 10) conditional_gwas

- **Purpose:** Run an approximate conditional association analysis from summary statistics and LD (COJO-style, exact OLS algebra).

- **Method:** The marginal correlation of each variant is r = z/sqrt(n - 2 + z^2). For every variant outside the conditioning set C, the standardised effect, standard error, t statistic and p of the variant are those of the joint OLS model {C + variant} computed exactly from the marginal correlations and the LD correlation matrix, with a t-test on n - |C| - 2 degrees of freedom. `r2_with_conditioning_set` is reported. The result is only as good as the LD matrix supplied (an in-sample LD reference is assumed).

- **API:** `dnakit.gwas.conditional_gwas(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. It needs `z` (or `beta` and `se`) and ideally `n`. `options`: `condition_on` (required; 1 to 20 variant ids), `ld` (`{'variant_ids': [...], 'matrix': [[...]]}` correlation matrix r) or `genotypes` (reference dosages), `n` (sample size; default the median `n` column).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `conditioned_on`, `n`, `records` (per other variant: `variant_id`, `marginal_z`, `marginal_p`, `conditional_beta_std`, `conditional_se_std`, `conditional_t`, `conditional_p`, `r2_with_conditioning_set`) and `method`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import conditional_gwas

summary = [
    {"variant_id": "v1", "beta": 0.30, "se": 0.05, "n": 1000},
    {"variant_id": "v2", "beta": 0.27, "se": 0.05, "n": 1000},
    {"variant_id": "v3", "beta": 0.05, "se": 0.05, "n": 1000},
]
ld = {"variant_ids": ["v1", "v2", "v3"],
      "matrix": [[1.0, 0.9, 0.1], [0.9, 1.0, 0.1], [0.1, 0.1, 1.0]]}
result = conditional_gwas(summary, options={"condition_on": ["v1"], "ld": ld})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "conditional_gwas",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 176,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "332f55a26240696b32413b21615aea7711c83f5056381b5d1ae715c26ff00db2"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "conditioned_on": [
      "v1"
    ],
    "n": 1000.0,
    "records": [
      {
        "variant_id": "v2",
        "marginal_z": 5.4,
        "marginal_p": 6.664089697085715e-08,
        "conditional_beta_std": 0.0029379654824893686,
        "conditional_se_std": 0.07138060961532508,
        "conditional_t": 0.04115915370185633,
        "conditional_p": 0.9671772573822426,
        "r2_with_conditioning_set": 0.81
      },
      {
        "variant_id": "v3",
        "marginal_z": 1.0,
        "marginal_p": 0.31731050786291415,
        "conditional_beta_std": 0.013110593385134293,
        "conditional_se_std": 0.031268103211206585,
        "conditional_t": 0.41929608894329784,
        "conditional_p": 0.6750900785228465,
        "r2_with_conditioning_set": 0.010000000000000009
      }
    ],
    "method": "exact joint OLS from marginal correlations and LD (COJO-style): standardised phenotype and genotypes, marginal r from z via r = z/sqrt(n-2+z^2), t-test with df = n - |C| - 2"
  }
}
```

## 11) define_association_loci

- **Purpose:** Define independent association loci as lead SNP +/- flank, merging overlapping windows.

- **Method:** Lead SNPs are selected as in `select_lead_snps` (`p_threshold`, `window`); each lead defines the interval [pos - flank, pos + flank] (start clamped at 1) and overlapping intervals on the same chromosome are merged into one locus. Each locus reports its leads, its strongest variant, the number of variants and the number of significant variants it contains.

- **API:** `dnakit.gwas.define_association_loci(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. It needs `chrom`, `pos` and `p`. `options`: `p_threshold` (5e-8), `window` (500000, used for lead selection), `flank` (500000).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `loci` (`chrom`, `start`, `end`, `leads`, `lead_variant`, `lead_p`, `n_variants`, `n_significant`), `n_loci` and `locus_definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import define_association_loci

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = define_association_loci(rows, options={"p_threshold": 5e-8, "flank": 100000})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "define_association_loci",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "loci": [
      {
        "chrom": "1",
        "start": 1,
        "end": 200000,
        "leads": [
          "rs1"
        ],
        "lead_variant": "rs1",
        "lead_p": 2e-09,
        "n_variants": 2,
        "n_significant": 1
      },
      {
        "chrom": "1",
        "start": 800000,
        "end": 1000000,
        "leads": [
          "rs3"
        ],
        "lead_variant": "rs3",
        "lead_p": 2e-09,
        "n_variants": 1,
        "n_significant": 1
      }
    ],
    "n_loci": 2,
    "locus_definition": "lead variant +/- flank bp (clamped at 1), overlapping loci on one chromosome merged"
  }
}
```

## 12) estimate_genetic_correlation

- **Purpose:** Estimate the SNP-based genetic correlation between two traits.

- **Method:** With individual data, cross-trait Haseman-Elston regression on the VanRaden GRM: with phenotypes scaled to unit variance, h2_1 and h2_2 are the slopes of y_i*y_j on GRM_ij over all pairs i<j (with intercept), the covariance is the slope of 0.5*(y1_i*y2_j + y1_j*y2_i), and rg = covariance/sqrt(h2_1*h2_2). With two summary-statistic tables, unweighted cross-trait LD-score regression: h2 per trait from the OLS of chi2 on the LD score, rho_g from the slope of z1*z2 on the LD score times M/sqrt(N1*N2) (the alleles of trait 2 are aligned to trait 1), rg = rho_g/sqrt(h1*h2). rg is `null` (with a warning) when a heritability estimate is not positive, and the estimate is not constrained to [-1, 1] (it is a ratio of noisy estimates and small samples can give a value such as 1.008): the value is kept unclamped and a warning is added whenever |rg| > 1. A warning is also added when an individual h2 estimate is below 0 or above 1, and the summary-statistics path always adds a warning that no standard errors are computed. No standard errors are produced. The GRM uses numpy when it is installed and a pure-Python path otherwise (refused above about 3e7 multiply-adds).

- **API:** `dnakit.gwas.estimate_genetic_correlation(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: either `{'genotypes': samples x variants dosages, 'phenotype_1': [...], 'phenotype_2': [...]}` or `{'trait1': table, 'trait2': table}` with `z` and LD scores (an `ldscore` column or `options['ld_scores']`; at least 10 variants). `options`: `ld_scores` ({variant_id: score}), `n1`, `n2`, `m`, `fix_intercept` (False).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. For individual data `result` holds `method` (`haseman_elston_cross_product`), `rg`, `genetic_covariance_std`, `h2_trait1`, `h2_trait2`, `n_samples` and `definition`; for summary statistics `method` (`ldsc_unweighted`), `rg`, `genetic_covariance`, `h2_trait1`, `h2_trait2`, `cross_intercept`, `n_variants` and `definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import estimate_genetic_correlation

genotypes = [
    [0, 1, 2, 0, 1, 2], [1, 0, 1, 2, 2, 0], [2, 2, 0, 1, 0, 1], [0, 1, 1, 0, 1, 2],
    [1, 2, 0, 2, 0, 1], [2, 0, 2, 1, 1, 0], [0, 2, 1, 1, 2, 2], [1, 1, 2, 0, 0, 1],
    [2, 1, 0, 2, 1, 0], [0, 0, 2, 2, 1, 1], [1, 2, 1, 0, 2, 0], [2, 0, 0, 1, 2, 2],
]
phenotype_1 = [1.2, 0.3, -0.8, 0.9, -0.1, 0.4, 1.0, 0.2, -1.1, 0.6, -0.4, -0.2]
phenotype_2 = [0.9, 0.5, -0.6, 0.7, -0.3, 0.1, 1.1, 0.0, -0.9, 0.8, -0.2, -0.5]
result = estimate_genetic_correlation(
    {"genotypes": genotypes, "phenotype_1": phenotype_1, "phenotype_2": phenotype_2}
)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "estimate_genetic_correlation",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 419,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "9d01b512cd38fea620a4c1050f752fb747559dbe156281ace6fa12cd20548bcd"
  },
  "warnings": ["rg = 1.008 is outside [-1, 1]: the ratio estimator is not bounded and small samples or noisy covariance can push it past 1; the value is reported unclamped and carries large uncertainty"],
  "status": "ok",
  "result": {
    "method": "haseman_elston_cross_product",
    "rg": 1.0079706026145199,
    "genetic_covariance_std": 0.8772411283478185,
    "h2_trait1": 0.9015109250540984,
    "h2_trait2": 0.8401778798532548,
    "n_samples": 12,
    "definition": "rg = slope(0.5*(y1_i*y2_j+y1_j*y2_i) ~ GRM_ij) / sqrt(h2_1*h2_2), phenotypes scaled to unit variance"
  }
}
```

## 13) estimate_snp_heritability

- **Purpose:** Estimate SNP heritability by Haseman-Elston regression (individual data) or unweighted LD-score regression (summary data).

- **Method:** `haseman_elston`: with the phenotype scaled to unit variance, h2 is the slope (with intercept) of y_i*y_j on the VanRaden GRM_ij over all pairs i<j; monomorphic markers are dropped. `ldsc`: OLS of chi2 = z^2 on the LD score l_j, E[chi2] = N*h2/M*l_j + intercept, so h2 = slope*M/N (the intercept is free unless `fix_intercept`, then 1); it needs at least 10 variants with z and an LD score. No standard errors or jackknife are produced and the estimates can be negative or above 1 in small samples. The value is kept unclamped, a warning is added whenever h2 < 0 or h2 > 1, and the summary-statistics (`ldsc`) path always adds a warning that no standard errors are computed. The GRM uses numpy when it is installed and a pure-Python path otherwise (refused above about 3e7 multiply-adds).

- **API:** `dnakit.gwas.estimate_snp_heritability(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'genotypes': samples x variants dosages, 'phenotype': [...]}` or a summary-statistics table with `z` and LD scores (`ld_score`/`ldscore` column or `options['ld_scores']`). `options`: `method` = `auto` (default; Haseman-Elston when `genotypes` is present, else `ldsc`), `haseman_elston` or `ldsc`; `ld_scores` ({variant_id: score}), `n` (default the median `n` column), `m` (default the number of variants used), `fix_intercept` (False).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. Haseman-Elston: `result` holds `method` (`haseman_elston_cross_product`), `h2`, `n_samples`, `n_markers` and `definition`. LD-score regression: `method` (`ldsc_unweighted`), `h2`, `intercept`, `slope`, `n`, `m`, `n_variants` and `definition`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import estimate_snp_heritability

genotypes = [
    [0, 1, 2, 0, 1, 2], [1, 0, 1, 2, 2, 0], [2, 2, 0, 1, 0, 1], [0, 1, 1, 0, 1, 2],
    [1, 2, 0, 2, 0, 1], [2, 0, 2, 1, 1, 0], [0, 2, 1, 1, 2, 2], [1, 1, 2, 0, 0, 1],
    [2, 1, 0, 2, 1, 0], [0, 0, 2, 2, 1, 1], [1, 2, 1, 0, 2, 0], [2, 0, 0, 1, 2, 2],
]
phenotype = [1.2, 0.3, -0.8, 0.9, -0.1, 0.4, 1.0, 0.2, -1.1, 0.6, -0.4, -0.2]
result = estimate_snp_heritability({"genotypes": genotypes, "phenotype": phenotype})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "estimate_snp_heritability",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 335,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f045666eb6204e967b3cc93c9d1e87ae7425f0fbf3754fb48f58fbbff9b878bf"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "method": "haseman_elston_cross_product",
    "h2": 0.9015109250540984,
    "n_samples": 12,
    "n_markers": 6,
    "definition": "slope of y_i*y_j on VanRaden GRM_ij over all pairs i<j (with intercept), y scaled to unit variance"
  }
}
```

## 14) evaluate_prs

- **Purpose:** Evaluate a polygenic score against a phenotype: R^2 and slope for quantitative traits, AUC and odds ratio for binary traits.

- **Method:** Quantitative trait: Pearson r, r^2 and the OLS slope of the phenotype on the standardised score (per SD) with a t-test p-value. Binary trait: AUC by the Mann-Whitney rank statistic (ties get mean ranks), and from a logistic regression on the standardised score the odds ratio per SD with a Wald 95% confidence interval and p-value and Nagelkerke R^2; if the logistic fit does not converge (perfect separation) only the AUC and a note are returned. Individuals with a missing value are dropped (at least 10 needed).

- **API:** `dnakit.gwas.evaluate_prs(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'score': [...], 'phenotype': [...]}`. `options`: `trait` = `auto` (default; binary when the phenotype is 0/1 with both classes), `quantitative` or `binary`.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. Quantitative: `result` holds `trait`, `n`, `r`, `r2`, `slope_per_sd`, `slope_se`, `p`. Binary: `trait`, `n`, `n_cases`, `n_controls`, `auc`, `odds_ratio_per_sd`, `or_ci95`, `log_or`, `se`, `p`, `nagelkerke_r2` (or a `note` if the logistic fit failed). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import evaluate_prs

score = [-1.2, -0.8, -0.5, -0.3, 0.0, 0.1, 0.3, 0.6, 0.9, 1.4, -1.0, -0.2, 0.2, 0.7, 1.1, 1.6]
phenotype = [0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 1]
result = evaluate_prs({"score": score, "phenotype": phenotype}, options={"trait": "binary"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "evaluate_prs",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 160,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "bf1d3cbdd5f0dd2548bfac81709c02e63edd6840e4f64c0e340776a98b2507e6"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "trait": "binary",
    "n": 16,
    "n_cases": 7,
    "n_controls": 9,
    "auc": 0.8888888888888888,
    "odds_ratio_per_sd": 10.592402008491819,
    "or_ci95": [
      1.1190821428442463,
      100.25982545332894
    ],
    "log_or": 2.360136952448344,
    "se": 1.1467701121929976,
    "p": 0.03958309732067795,
    "nagelkerke_r2": 0.5816456123759354
  }
}
```

## 15) export

- **Purpose:** Serialise GWAS summary statistics to another file format (tsv, csv, json, jsonl, bed, ldsc, metal).

- **Method:** The table is standardised and written in the chosen layout: `tsv`/`csv` (all non-empty canonical columns, or `columns`), `json`, `jsonl`, `bed` (chr-prefixed chromosome, 0-based start, end, variant id and -log10 p as the score), `ldsc` (the munged `SNP A1 A2 N Z` layout) or `metal` (`MarkerName Allele1 Allele2 Freq1 Effect StdErr P-value`). The text is always returned; a file is written only when `path` is given (`.gz` is gzip-compressed) and an existing file is not overwritten unless `overwrite=True`.

- **API:** `dnakit.gwas.export(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. `options`: `format` = `tsv` (default), `csv`, `json`, `jsonl`, `bed`, `ldsc` or `metal`; `columns` (subset for tsv/csv); `path` (or `output`); `overwrite` (False). `bed` needs `chrom` and `pos`; `ldsc` needs `variant_id`, both alleles and `z`; `metal` needs `variant_id`, both alleles, `beta` and `se`.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `format`, `n_records`, `text` (the serialised table) and `path` (null unless written). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import export

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = export(rows[:3], options={"format": "metal"})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "export",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 374,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "879d97b1e5ae6c404c3f251cf85e86d9b0381c2aa31d5343b9285bc0814d8299"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "metal",
    "n_records": 3,
    "text": "MarkerName\tAllele1\tAllele2\tFreq1\tEffect\tStdErr\tP-value\nrs1\tA\tG\t0.3\t0.12\t0.02\t2e-09\nrs2\tC\tT\t0.25\t0.09\t0.02\t6.8e-06\nrs3\tG\tA\t0.4\t-0.15\t0.025\t2e-09\n",
    "path": null
  }
}
```

## 16) filter

- **Purpose:** Filter GWAS summary-statistics rows by p-value, MAF, sample size, INFO, chromosome, region or variant type.

- **Method:** Each row is tested against the active filters in a fixed order (`require_p`, `max_p`, `min_maf` with MAF = min(eaf, 1-eaf), `min_n`, `min_info`, `chromosomes`, `region`, `snp_only`, `remove_palindromic`); a row with a missing value for an active filter is removed, and a removed row is counted under the first reason that applies.

- **API:** `dnakit.gwas.filter(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. `options`: `max_p`, `min_maf`, `min_n`, `min_info`, `chromosomes` (list), `region` (`chrom:start-end`), `snp_only` (False), `remove_palindromic` (False), `require_p` (False). Filters that are not given are not applied.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `records` (the kept, standardised rows), `n_input`, `n_kept` and `removed` ({reason: count}, with reasons such as `p_above_threshold`, `maf_below_threshold`, `n_below_threshold`, `info_below_threshold`, `chromosome_excluded`, `outside_region`, `not_a_snp`, `palindromic`, `missing_p`). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import filter

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = filter(rows, options={"max_p": 0.05, "min_maf": 0.2, "chromosomes": ["1", "2"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "gwas",
  "operation": "filter",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "records": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.12,
        "se": 0.02,
        "p": 2e-09,
        "eaf": 0.3,
        "n": 5000,
        "z": 6.0
      },
      {
        "variant_id": "rs2",
        "chrom": "1",
        "pos": 150000,
        "effect_allele": "C",
        "other_allele": "T",
        "beta": 0.09,
        "se": 0.02,
        "p": 6.8e-06,
        "eaf": 0.25,
        "n": 5000,
        "z": 4.5
      },
      "..."
    ],
    "n_input": 6,
    "n_kept": 4,
    "removed": {
      "p_above_threshold": 2
    }
  }
}
```

## 17) fixed_effect_meta_analysis

- **Purpose:** Combine per-variant effects across studies with inverse-variance-weighted fixed-effect meta-analysis.

- **Method:** For each variant (matched by `variant_id`, or `chrom:pos`), weights w = 1/se^2, pooled beta = sum(w*beta)/sum(w), se = 1/sqrt(sum(w)), z = beta/se with a two-sided normal p-value, Cochran's Q = sum(w*(beta - pooled)^2) with a chi-square p-value on k-1 df, and I^2 = max(0, (Q - df)/Q). The alleles of every study are aligned to the first study that reports the variant (the effect is negated when the alleles are swapped; studies with mismatching alleles are skipped, strand flips are not resolved). Variants seen in fewer than `min_studies` studies are omitted.

- **API:** `dnakit.gwas.fixed_effect_meta_analysis(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a dict {study name: table} or a list of tables (each with `variant_id` or `chrom` + `pos`, `beta`, `se`, alleles), or one table with a `study` column. `options`: `min_studies` (2), `study_column` (`study`). At least two studies are needed.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `records` (per variant: `variant_id`, `chrom`, `pos`, `effect_allele`, `other_allele`, `n_studies`, `beta`, `se`, `z`, `p`, `q`, `q_p`, `i2`), `n_variants` and `model`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import fixed_effect_meta_analysis

study_a = [{"variant_id": "rs1", "effect_allele": "A", "other_allele": "G", "beta": 0.10, "se": 0.03},
           {"variant_id": "rs2", "effect_allele": "C", "other_allele": "T", "beta": 0.02, "se": 0.03}]
study_b = [{"variant_id": "rs1", "effect_allele": "A", "other_allele": "G", "beta": 0.14, "se": 0.04},
           {"variant_id": "rs2", "effect_allele": "T", "other_allele": "C", "beta": -0.05, "se": 0.04}]
result = fixed_effect_meta_analysis({"A": study_a, "B": study_b})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "fixed_effect_meta_analysis",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 382,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "3113603d7296e5791a4642b96f6fd8c8a8dd6407e46ab7d703c2c4f3b959d454"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "records": [
      {
        "variant_id": "rs1",
        "chrom": null,
        "pos": null,
        "effect_allele": "A",
        "other_allele": "G",
        "n_studies": 2,
        "beta": 0.11440000000000002,
        "se": 0.024,
        "z": 4.7666666666666675,
        "p": 1.872988003352597e-06,
        "q": 0.6400000000000003,
        "q_p": 0.4237107971667936,
        "i2": 0.0
      },
      {
        "variant_id": "rs2",
        "chrom": null,
        "pos": null,
        "effect_allele": "C",
        "other_allele": "T",
        "n_studies": 2,
        "beta": 0.0308,
        "se": 0.024,
        "z": 1.2833333333333334,
        "p": 0.19937531925904434,
        "q": 0.36,
        "q_p": 0.5485062355001472,
        "i2": 0.0
      }
    ],
    "n_variants": 2,
    "model": "fixed effect, inverse-variance weights w=1/se^2; Cochran's Q and I^2 reported"
  }
}
```

## 18) harmonize

- **Purpose:** Harmonise the effect alleles of a summary-statistics table to a reference table (allele swap, strand flip and palindrome handling).

- **Method:** Variants are matched to the reference by `chrom:pos` (or by id; `match_on`='auto' uses position when every row has it). If the effect allele is the reference's other allele the alleles are swapped and beta and z are negated and eaf becomes 1-eaf; alleles on the opposite strand are complemented (`strand_flip`, `strand_flip_swap`). Non-palindromic mismatches are dropped. Palindromic (A/T, C/G) SNPs are ambiguous: with `palindromes='drop_ambiguous'` they are dropped when the target or the reference frequency is within `palindrome_eaf_tolerance` of 0.5 and otherwise resolved by comparing the frequencies; with `assume_forward_strand` they are kept without strand inference. Every kept row carries `harmonize_action`.

- **API:** `dnakit.gwas.harmonize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'target': table, 'reference': table}`, or a table with `options['reference']`; both need effect and other alleles. `options`: `match_on` = `auto` (default), `id` or `position`; `palindromes` = `drop_ambiguous` (default) or `assume_forward_strand`; `palindrome_eaf_tolerance` (0.08).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `records` (harmonised rows with `harmonize_action`: `none`, `swap`, `strand_flip`, `strand_flip_swap`), `n_target`, `n_harmonized`, `actions` (counts, including the drop reasons `not_in_reference`, `missing_alleles`, `allele_mismatch`, `ambiguous_palindrome`) and `matched_on`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import harmonize

target = [
    {"variant_id": "rs1", "chrom": "1", "pos": 100, "effect_allele": "A", "other_allele": "G", "beta": 0.10, "eaf": 0.30},
    {"variant_id": "rs2", "chrom": "1", "pos": 200, "effect_allele": "T", "other_allele": "C", "beta": 0.20, "eaf": 0.20},
    {"variant_id": "rs3", "chrom": "1", "pos": 300, "effect_allele": "A", "other_allele": "T", "beta": 0.05, "eaf": 0.49},
]
reference = [
    {"variant_id": "rs1", "chrom": "1", "pos": 100, "effect_allele": "A", "other_allele": "G"},
    {"variant_id": "rs2", "chrom": "1", "pos": 200, "effect_allele": "C", "other_allele": "T"},
    {"variant_id": "rs3", "chrom": "1", "pos": 300, "effect_allele": "A", "other_allele": "T"},
]
result = harmonize({"target": target, "reference": reference})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "harmonize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 656,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "cf206a3f2cac9b758aa78197bc7656674394c05b3e6cc4872f65ad7220010816"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "records": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.1,
        "se": null,
        "p": null,
        "eaf": 0.3,
        "n": null,
        "z": null,
        "harmonize_action": "none"
      },
      {
        "variant_id": "rs2",
        "chrom": "1",
        "pos": 200,
        "effect_allele": "C",
        "other_allele": "T",
        "beta": -0.2,
        "se": null,
        "p": null,
        "eaf": 0.8,
        "n": null,
        "z": null,
        "harmonize_action": "swap"
      }
    ],
    "n_target": 3,
    "n_harmonized": 2,
    "actions": {
      "none": 1,
      "swap": 1,
      "ambiguous_palindrome": 1
    },
    "matched_on": "position"
  }
}
```

## 19) joint_gwas_analysis

- **Purpose:** Run a joint multi-SNP association model from summary statistics and LD (COJO-joint style).

- **Method:** Computes the exact multiple regression of the standardised phenotype on the selected standardised genotypes from the marginal correlations r = z/sqrt(n - 2 + z^2) and the LD correlation matrix R: beta = R^-1 r, model R^2 = r'beta, se^2 = (1 - R^2)/(n - k - 1) * diag(R^-1) and a t-test on n - k - 1 degrees of freedom. A singular LD matrix (perfect LD) is rejected. The result is only as good as the LD matrix supplied (an in-sample LD reference is assumed).

- **API:** `dnakit.gwas.joint_gwas_analysis(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. It needs `z` (or `beta` and `se`) and ideally `n`. `options`: `ld` (`{'variant_ids': [...], 'matrix': [[...]]}` correlation matrix r) or `genotypes` (reference dosages), `snps` (variant ids to model jointly, at most 200; default all in the LD matrix that are in the table), `n` (default the median `n` column).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `n`, `records` (per variant: `variant_id`, `marginal_z`, `joint_beta_std`, `joint_se_std`, `joint_t`, `joint_p`), `model_r2`, `df` and `method`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import joint_gwas_analysis

summary = [
    {"variant_id": "v1", "beta": 0.30, "se": 0.05, "n": 1000},
    {"variant_id": "v2", "beta": 0.27, "se": 0.05, "n": 1000},
    {"variant_id": "v3", "beta": 0.05, "se": 0.05, "n": 1000},
]
ld = {"variant_ids": ["v1", "v2", "v3"],
      "matrix": [[1.0, 0.9, 0.1], [0.9, 1.0, 0.1], [0.1, 0.1, 1.0]]}
result = joint_gwas_analysis(summary, options={"ld": ld})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "joint_gwas_analysis",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 176,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "332f55a26240696b32413b21615aea7711c83f5056381b5d1ae715c26ff00db2"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "n": 1000.0,
    "records": [
      {
        "variant_id": "v1",
        "marginal_z": 5.999999999999999,
        "joint_beta_std": 0.1832581196160083,
        "joint_se_std": 0.07142915304024819,
        "joint_t": 2.565592784122021,
        "joint_p": 0.01044535198931147
      },
      {
        "variant_id": "v2",
        "marginal_z": 5.4,
        "joint_beta_std": 0.002249129961401298,
        "joint_se_std": 0.07142915304024819,
        "joint_t": 0.03148756307013721,
        "joint_p": 0.9748870188760195
      },
      {
        "variant_id": "v3",
        "marginal_z": 1.0,
        "joint_beta_std": 0.013087874900675695,
        "joint_se_std": 0.03129209958977394,
        "joint_t": 0.4182485378818342,
        "joint_p": 0.6758555722527371
      }
    ],
    "model_r2": 0.03498737698481886,
    "df": 996.0,
    "method": "exact multiple regression of the standardised phenotype on the selected standardised genotypes, computed from marginal correlations and LD (beta = R^-1 r, se^2 = (1-R2)/(n-k-1) * diag(R^-1))"
  }
}
```

## 20) merge_nearby_significant_hits

- **Purpose:** Merge nearby genome-wide-significant variants into hit regions by single-linkage on base-pair distance.

- **Method:** The variants with p <= `p_threshold` are sorted by chromosome and position; consecutive significant variants on one chromosome at most `distance` bp apart share a region. Each region reports its start and end positions, the number of hits and its lead (smallest-p) variant.

- **API:** `dnakit.gwas.merge_nearby_significant_hits(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. It needs `chrom`, `pos` and `p`. `options`: `p_threshold` (5e-8), `distance` (500000).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `groups` (`chrom`, `start`, `end`, `n_hits`, `lead_variant`, `lead_p`), `n_groups`, `n_hits`, `distance` and `rule`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import merge_nearby_significant_hits

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = merge_nearby_significant_hits(rows, options={"p_threshold": 1e-5, "distance": 100000})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "merge_nearby_significant_hits",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "groups": [
      {
        "chrom": "1",
        "start": 100000,
        "end": 150000,
        "n_hits": 2,
        "lead_variant": "rs1",
        "lead_p": 2e-09
      },
      {
        "chrom": "1",
        "start": 900000,
        "end": 900000,
        "n_hits": 1,
        "lead_variant": "rs3",
        "lead_p": 2e-09
      }
    ],
    "n_groups": 2,
    "n_hits": 3,
    "distance": 100000,
    "rule": "single-linkage: consecutive significant variants closer than or equal to `distance` bp are merged"
  }
}
```

## 21) plot

- **Purpose:** Draw a Manhattan or QQ plot of GWAS p-values as an SVG string.

- **Method:** The points come from `compute_manhattan_data` or `compute_qq_data`. Manhattan: alternating chromosome colours, dashed lines at the genome-wide (5e-8) and suggestive (1e-5) thresholds and chromosome labels at the label centres. QQ: observed against expected -log10 p with the y = x line. Pure standard-library SVG; nothing is written unless `path` is given.

- **API:** `dnakit.gwas.plot(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. `options`: `kind` = `manhattan` (default; needs `chrom`, `pos`, `p`) or `qq` (needs `p`); `width` (900), `height` (450), `max_points` (20000), `path` (an existing file is not overwritten unless `overwrite=True`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `kind`, `svg` (the SVG text), `n_points_drawn` and `path` (null unless written). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import plot

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = plot(rows, options={"kind": "manhattan"})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "gwas",
  "operation": "plot",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "kind": "manhattan",
    "svg": "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"900\" height=\"450\" viewBox=\"0 0 900 450\"><title>GWAS manhattan plot</title...",
    "n_points_drawn": 6,
    "path": null
  }
}
```

## 22) random_effect_meta_analysis

- **Purpose:** Combine per-variant effects across studies with DerSimonian-Laird random-effects meta-analysis.

- **Method:** Input handling and allele alignment are as in `fixed_effect_meta_analysis`. With w = 1/se^2, S1 = sum(w), S2 = sum(w^2) and Cochran's Q on k-1 df, tau^2 = max(0, (Q - df)/(S1 - S2/S1)); the pooled effect uses weights 1/(se^2 + tau^2), se = 1/sqrt(sum of these weights) and a normal-theory p-value. Q, its p-value and I^2 are those of the fixed-effect model. With few studies tau^2 is poorly estimated.

- **API:** `dnakit.gwas.random_effect_meta_analysis(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** As `fixed_effect_meta_analysis`: a dict {study name: table}, a list of tables, or one table with a `study` column; `options`: `min_studies` (2), `study_column` (`study`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `records` (per variant: `variant_id`, `chrom`, `pos`, `effect_allele`, `other_allele`, `n_studies`, `beta`, `se`, `z`, `p`, `tau2`, `q`, `q_p`, `i2`), `n_variants` and `model`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import random_effect_meta_analysis

studies = {
    "A": [{"variant_id": "rs1", "effect_allele": "A", "other_allele": "G", "beta": 0.10, "se": 0.03}],
    "B": [{"variant_id": "rs1", "effect_allele": "A", "other_allele": "G", "beta": 0.30, "se": 0.04}],
    "C": [{"variant_id": "rs1", "effect_allele": "A", "other_allele": "G", "beta": 0.05, "se": 0.05}],
}
result = random_effect_meta_analysis(studies)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "random_effect_meta_analysis",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 295,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "85abbef3f8debc7c12e8a2bc7be9bb7a242a534e0d4eee0f114d34a321caef77"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "records": [
      {
        "variant_id": "rs1",
        "chrom": null,
        "pos": null,
        "effect_allele": "A",
        "other_allele": "G",
        "n_studies": 3,
        "beta": 0.15101036270181903,
        "se": 0.07331952361816618,
        "z": 2.05962007456911,
        "p": 0.03943487558282798,
        "tau2": 0.014487000000000002,
        "q": 20.838751625487642,
        "q_p": 2.9848504713101197e-05,
        "i2": 0.9040249609984399
      }
    ],
    "n_variants": 1,
    "model": "DerSimonian-Laird random effects: tau^2 = max(0, (Q-df)/(S1 - S2/S1)), weights 1/(se^2+tau^2), normal-theory p"
  }
}
```

## 23) read

- **Purpose:** Read a GWAS summary-statistics table (path, .gz file, delimited text, records or columns) into records.

- **Method:** The delimiter (tab, comma, semicolon or whitespace) is detected, `NA`/empty cells become null and numeric text becomes numbers. By default the original column names are kept; with `standardize=True` the columns are renamed to the canonical names (see `standardize`). Tables of more than `max_rows` rows are refused.

- **API:** `dnakit.gwas.read(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: a file path (`.gz` allowed), delimited text, a list of records or a dict of columns. `options`: `delimiter`, `max_rows` (500000), `standardize` (False).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `columns`, `n_records` and `records`; with `standardize=True` also `column_map`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import read

text = "SNP\tCHR\tBP\tA1\tA2\tBETA\tSE\tP\nrs1\t1\t100000\tA\tG\t0.12\t0.02\t2e-9\nrs2\t1\t150000\tC\tT\t0.09\t0.02\tNA\n"
result = read(text)
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "gwas",
  "operation": "read",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 3,
    "text_length": 89,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "ce0498fb00c9e38118484def9a109a673112328311015fc7e7c629675152e6dd"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "columns": [
      "SNP",
      "CHR",
      "BP",
      "..."
    ],
    "n_records": 2,
    "records": [
      {
        "SNP": "rs1",
        "CHR": 1,
        "BP": 100000,
        "A1": "A",
        "A2": "G",
        "BETA": 0.12,
        "SE": 0.02,
        "P": 2e-09
      },
      {
        "SNP": "rs2",
        "CHR": 1,
        "BP": 150000,
        "A1": "C",
        "A2": "T",
        "BETA": 0.09,
        "SE": 0.02,
        "P": null
      }
    ]
  }
}
```

## 24) run

- **Purpose:** Run a native single-variant association scan (linear or logistic regression) on dosage genotypes.

- **Method:** Per variant, with an intercept and the covariates, a quantitative phenotype is tested by OLS (t-test on the dosage coefficient, df = n - covariates - 2) and a 0/1 phenotype by logistic regression (Wald z-test); `model='auto'` picks logistic for a 0/1 phenotype with both classes. Samples with a missing genotype or phenotype are dropped per variant. Monomorphic variants or too few samples get `status` `monomorphic_or_too_few_samples`, variants below `maf_min` `maf_below_threshold`, a failed fit `model_failed`. It is a plain regression scan: no relatedness or population-structure correction beyond the covariates supplied. External programs (PLINK, REGENIE, GEMMA, GCTA) are not driven: with `backend` set to one of them the call returns `not_implemented`, naming the missing program, or, when the program is installed, stating that no command-line adapter exists in this release.

- **API:** `dnakit.gwas.run(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: `{'genotypes': samples x variants dosages in [0, 2] (null = missing), 'phenotype': [...], 'covariates': samples x k (optional), 'variant_ids', 'chrom', 'pos', 'alt', 'ref', 'sample_ids'}` (at least 5 samples with a phenotype). `options`: `model` = `auto` (default), `linear` or `logistic`; `maf_min` (0); `backend` = `native` (default) or `plink`, `plink2`, `plink1.9`, `regenie`, `gemma`, `gcta` (these report `not_implemented`).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `model` (`linear` or `logistic`), `n_samples`, `n_variants`, `records` (per variant: `variant_id`, `chrom`, `pos`, `effect_allele`, `other_allele`, `beta`, `se`, `z`, `p`, `eaf`, `n`, `status`) and `test`. With an external `backend` the envelope has `status: "not_implemented"`, `required_dependency` (the program name, e.g. `plink2`) and `message`; an unknown backend name raises an error. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import run

genotypes = [[0, 1, 2], [1, 1, 0], [2, 0, 1], [0, 2, 1], [1, 0, 2], [2, 1, 0], [0, 0, 1], [1, 2, 2]]
phenotype = [0.1, 0.9, 1.8, 0.2, 1.1, 2.0, 0.0, 1.2]
result = run({"genotypes": genotypes, "phenotype": phenotype,
              "variant_ids": ["v1", "v2", "v3"], "chrom": ["1", "1", "1"], "pos": [100, 200, 300]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "run",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 1,
    "text_length": 243,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "7e3665d19858b0c6bbf219dda78749a79278644f55aaaebfdafbd34aedfbebcf"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "model": "linear",
    "n_samples": 8,
    "n_variants": 3,
    "records": [
      {
        "variant_id": "v1",
        "chrom": "1",
        "pos": 100,
        "effect_allele": null,
        "other_allele": null,
        "beta": 0.905128205128205,
        "se": 0.05695155982549665,
        "z": 15.892948461843325,
        "p": 3.938267710360089e-06,
        "eaf": 0.4375,
        "n": 8,
        "status": "ok"
      },
      {
        "variant_id": "v2",
        "chrom": "1",
        "pos": 200,
        "effect_allele": null,
        "other_allele": null,
        "beta": -0.12051282051282054,
        "se": 0.3706288819493145,
        "z": -0.32515766142936836,
        "p": 0.7561008935006274,
        "eaf": 0.4375,
        "n": 8,
        "status": "ok"
      },
      {
        "variant_id": "v3",
        "chrom": "1",
        "pos": 300,
        "effect_allele": null,
        "other_allele": null,
        "beta": -0.2897435897435897,
        "se": 0.35467510976489625,
        "z": -0.8169267641466362,
        "p": 0.4451878633327042,
        "eaf": 0.5625,
        "n": 8,
        "status": "ok"
      }
    ],
    "test": "OLS t-test on the dosage coefficient (df = n - covariates - 2)"
  }
}
```

## 25) select_lead_snps

- **Purpose:** Select independent lead SNPs by greedy distance-based pruning of genome-wide-significant variants.

- **Method:** Variants with p <= `p_threshold` are taken in order of ascending p; each becomes a lead and removes the other candidates within `window` bp on the same chromosome (their number is `n_nearby_significant`). No LD is used. Leads are returned in genomic order; `max_leads` stops the selection early.

- **API:** `dnakit.gwas.select_lead_snps(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. It needs `chrom`, `pos` and `p`. `options`: `p_threshold` (5e-8), `window` (500000), `max_leads`.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `leads` (`variant_id`, `chrom`, `pos`, alleles, `beta`, `se`, `p`, `eaf`, `neglog10p`, `n_nearby_significant`), `n_leads`, `p_threshold` and `window`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import select_lead_snps

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = select_lead_snps(rows, options={"p_threshold": 5e-8, "window": 500000})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "select_lead_snps",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 6,
    "text_length": 742,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "038644c8ea1ee6bddafba22cdada5022e0f55566e4c25e56c063727d08bdb0ff"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "leads": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.12,
        "se": 0.02,
        "p": 2e-09,
        "eaf": 0.3,
        "neglog10p": 8.698970004336019,
        "n_nearby_significant": 0
      },
      {
        "variant_id": "rs3",
        "chrom": "1",
        "pos": 900000,
        "effect_allele": "G",
        "other_allele": "A",
        "beta": -0.15,
        "se": 0.025,
        "p": 2e-09,
        "eaf": 0.4,
        "neglog10p": 8.698970004336019,
        "n_nearby_significant": 0
      }
    ],
    "n_leads": 2,
    "p_threshold": 5e-08,
    "window": 500000
  }
}
```

## 26) simulate_gwas_data

- **Purpose:** Simulate genotypes, a quantitative phenotype and the resulting association summary statistics (deterministic by seed).

- **Method:** Allele frequencies f_j are uniform on [maf_min, maf_max] and dosages are Binomial(2, f_j) (unlinked variants). `n_causal` randomly chosen variants get effects b_j ~ N(0, heritability/n_causal) on the standardised dosage; phenotype = sum_j b_j * standardised dosage_j + N(0, 1 - heritability). Variants are spread over `n_chromosomes` chromosomes at 10 kb spacing with alt allele `A` and ref `G`. The summary statistics are the linear-regression scan of `run` on the simulated data. Python `random`, deterministic given `seed`; there is no LD structure.

- **API:** `dnakit.gwas.simulate_gwas_data(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** `input`: None or a mapping of the options below (it overrides `options`). `options`: `n_samples` (500), `n_variants` (200), `n_causal` (5), `heritability` (0.5), `seed` (0), `maf_min` (0.05), `maf_max` (0.5), `n_chromosomes` (22).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `genotypes` (a genotype mapping usable by `run` and the other functions: `genotypes`, `sample_ids`, `variant_ids`, `chrom`, `pos`, `alt`, `ref`), `phenotype`, `true_effects` ({variant_id: effect}), `allele_frequencies`, `summary_statistics` (the `run` records), `seed` and `model`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import simulate_gwas_data

result = simulate_gwas_data(None, options={"n_samples": 50, "n_variants": 10, "n_causal": 2, "seed": 1})
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "gwas",
  "operation": "simulate_gwas_data",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 0,
    "text_length": 0,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "genotypes": {
      "genotypes": [
        [
          0,
          1,
          0,
          "..."
        ],
        [
          0,
          1,
          2,
          "..."
        ],
        "..."
      ],
      "sample_ids": [
        "s1",
        "s2",
        "s3",
        "..."
      ],
      "variant_ids": [
        "snp1",
        "snp2",
        "snp3",
        "..."
      ],
      "chrom": [
        "1",
        "2",
        "3",
        "..."
      ],
      "pos": [
        10000,
        10000,
        10000,
        "..."
      ],
      "alt": [
        "A",
        "A",
        "A",
        "..."
      ],
      "ref": [
        "G",
        "G",
        "G",
        "..."
      ]
    },
    "phenotype": [
      -0.04701902898502952,
      0.9976611520850311,
      -0.618585574379673,
      "..."
    ],
    "true_effects": {
      "snp5": -0.132523395903016,
      "snp8": -0.3711891596172627
    },
    "allele_frequencies": {
      "snp1": 0.11046390985058055,
      "snp2": 0.4313451816217547,
      "snp3": 0.3936985785394763,
      "snp4": 0.16478106158273975,
      "snp5": 0.2729457891913734,
      "snp6": 0.2522709791549322,
      "snp7": 0.3432168377252433,
      "snp8": 0.40492550801098093,
      "snp9": 0.0922368140484057,
      "snp10": 0.06275636443490284
    },
    "summary_statistics": [
      {
        "variant_id": "snp1",
        "chrom": "1",
        "pos": 10000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": -0.06099515839356878,
        "se": 0.18813268750275483,
        "z": -0.3242135069838708,
        "p": 0.7471850214511193,
        "eaf": 0.14,
        "n": 50,
        "status": "ok"
      },
      {
        "variant_id": "snp2",
        "chrom": "2",
        "pos": 10000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": -0.0013527033652641059,
        "se": 0.14615643899646016,
        "z": -0.00925517462352013,
        "p": 0.9926539046237707,
        "eaf": 0.41,
        "n": 50,
        "status": "ok"
      },
      "..."
    ],
    "seed": 1,
    "model": "dosage ~ Binomial(2, f_j); phenotype = sum_j b_j * standardised dosage_j + N(0, 1 - h2), b_j ~ N(0, h2 / n_causal)"
  }
}
```

## 27) standardize

- **Purpose:** Standardise GWAS summary-statistics columns to the canonical names and derive missing statistics.

- **Method:** Recognised aliases are renamed to `variant_id`, `chrom` (no `chr` prefix; 23 -> X, 24 -> Y, 26 -> MT), `pos`, `effect_allele` and `other_allele` (upper case), `beta`, `se`, `p`, `eaf`, `n` and `z`; `maf`, `info` and `ld_score` columns are kept when present. Missing values are derived: beta = ln(OR) from an odds ratio, z = beta/se, p from z (two-sided normal), se from beta and p, and `variant_id` = `chrom:pos:other:effect` when absent. Unrecognised columns are dropped with a warning.

- **API:** `dnakit.gwas.standardize(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. `options`: `max_rows` (500000), `delimiter`.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `records` (rows with the canonical columns), `column_map` ({original name: canonical name}), `columns` and `n_records`; `warnings` lists the ignored columns. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import standardize

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = standardize(rows[:2])
print(json.dumps(result, indent=2))
```

- **Example output:**

Long arrays and strings are shortened with `...`.

```json
{
  "module": "gwas",
  "operation": "standardize",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 249,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a530a3275efbe926d438abe1d12407274b4c4cbfd0513d44cad410c43b350a28"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "records": [
      {
        "variant_id": "rs1",
        "chrom": "1",
        "pos": 100000,
        "effect_allele": "A",
        "other_allele": "G",
        "beta": 0.12,
        "se": 0.02,
        "p": 2e-09,
        "eaf": 0.3,
        "n": 5000,
        "z": 6.0
      },
      {
        "variant_id": "rs2",
        "chrom": "1",
        "pos": 150000,
        "effect_allele": "C",
        "other_allele": "T",
        "beta": 0.09,
        "se": 0.02,
        "p": 6.8e-06,
        "eaf": 0.25,
        "n": 5000,
        "z": 4.5
      }
    ],
    "column_map": {
      "SNP": "variant_id",
      "CHR": "chrom",
      "BP": "pos",
      "A1": "effect_allele",
      "A2": "other_allele",
      "BETA": "beta",
      "SE": "se",
      "P": "p",
      "EAF": "eaf",
      "N": "n"
    },
    "columns": [
      "variant_id",
      "chrom",
      "pos",
      "..."
    ],
    "n_records": 2
  }
}
```

## 28) validate

- **Purpose:** Validate a GWAS summary-statistics table: required columns, value ranges, alleles and duplicate variants.

- **Method:** After standardisation it reports errors (no variant identifier or `chrom` + `pos`; no `p`, `z` or `beta` with `se`; p outside [0, 1]; se <= 0; eaf outside [0, 1]; n <= 0; non-finite beta/se/z; identical effect and other alleles) and warnings (alleles other than A/C/G/T strings, duplicated variant ids, missing allele or eaf columns). Content problems give `valid: false` rather than raising; an unreadable table raises.

- **API:** `dnakit.gwas.validate(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. `options`: `max_rows` (500000), `delimiter`.

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `valid`, `n_records`, `errors` (list of messages), `warnings` (list of messages) and `column_map`. Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import validate

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100, "A1": "A", "A2": "G", "BETA": 0.1, "SE": 0.02, "P": 1.5},
    {"SNP": "rs1", "CHR": 1, "BP": 200, "A1": "C", "A2": "C", "BETA": 0.2, "SE": -0.02, "P": 0.01},
]
result = validate(rows)
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "validate",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 190,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "f4830ef0886039ab4da0ddc10686680f90e3be2a0f403baedda6e82657802523"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "valid": false,
    "n_records": 2,
    "errors": [
      "1 p-values outside [0, 1]",
      "1 standard errors that are not > 0",
      "1 variants have identical effect and other alleles"
    ],
    "warnings": [
      "1 duplicated variant identifiers",
      "no effect-allele frequency column"
    ],
    "column_map": {
      "SNP": "variant_id",
      "CHR": "chrom",
      "BP": "pos",
      "A1": "effect_allele",
      "A2": "other_allele",
      "BETA": "beta",
      "SE": "se",
      "P": "p"
    }
  }
}
```

## 29) write

- **Purpose:** Write GWAS summary statistics as a delimited table (tsv or csv) with canonical column names.

- **Method:** The table is standardised and serialised with the non-empty canonical columns, or with `columns`. The text is always returned; a file is written only when `path` is given (`.gz` is gzip-compressed) and an existing file is not overwritten unless `overwrite=True`. The format follows `format`, else the path suffix (`.csv` gives csv, anything else tsv).

- **API:** `dnakit.gwas.write(input[required], options[optional])` (formal DNAKit namespace; Agent/MCP exposure depends on the JSON call plan).

- **Input:** Summary statistics: a list of dict records, a dict of equal-length columns, or delimited text / a path (`.gz` allowed). Column names are matched case- and punctuation-insensitively to `variant_id` (SNP, rsid), `chrom` (CHR), `pos` (BP, POS), `effect_allele` (A1, EA, ALT), `other_allele` (A2, NEA, REF), `beta`, `se`, `p` (P, PVAL), `eaf` (EAF, FREQ), `n`, `z` and `OR` (converted to beta = ln OR); missing `z`, `p`, `se` and `variant_id` are derived where possible. `options`: `path` (or `output`), `format` = `tsv` or `csv`, `columns` (subset, canonical names), `overwrite` (False).

- **Output:** Returns the DNAKit envelope: `status` (`ok`), `module`, `operation`, `provenance` (`implementation: native`), `input_summary`, `warnings` and `result`. `result` holds `format`, `columns`, `n_records`, `text` and `path` (null unless a file was written). Malformed input or invalid options raise a DNAKit exception.

- **Example code:**

```python
import json
from dnakit.gwas import write

rows = [
    {"SNP": "rs1", "CHR": 1, "BP": 100000, "A1": "A", "A2": "G", "BETA": 0.12, "SE": 0.02, "P": 2.0e-9, "EAF": 0.30, "N": 5000},
    {"SNP": "rs2", "CHR": 1, "BP": 150000, "A1": "C", "A2": "T", "BETA": 0.09, "SE": 0.02, "P": 6.8e-6, "EAF": 0.25, "N": 5000},
    {"SNP": "rs3", "CHR": 1, "BP": 900000, "A1": "G", "A2": "A", "BETA": -0.15, "SE": 0.025, "P": 2.0e-9, "EAF": 0.40, "N": 5000},
    {"SNP": "rs4", "CHR": 2, "BP": 50000, "A1": "T", "A2": "C", "BETA": 0.01, "SE": 0.02, "P": 0.62, "EAF": 0.10, "N": 5000},
    {"SNP": "rs5", "CHR": 2, "BP": 70000, "A1": "A", "A2": "C", "BETA": 0.04, "SE": 0.02, "P": 0.0455, "EAF": 0.50, "N": 5000},
    {"SNP": "rs6", "CHR": 2, "BP": 90000, "A1": "G", "A2": "T", "BETA": -0.02, "SE": 0.02, "P": 0.317, "EAF": 0.35, "N": 5000},
]
result = write(rows[:2], options={"format": "tsv", "columns": ["variant_id", "chrom", "pos", "beta", "p"]})
print(json.dumps(result, indent=2))
```

- **Example output:**

```json
{
  "module": "gwas",
  "operation": "write",
  "provenance": {
    "implementation": "native",
    "backend": "local"
  },
  "input_summary": {
    "record_count": 2,
    "text_length": 249,
    "sequence_bases": 0,
    "lengths": [],
    "sha256": "a530a3275efbe926d438abe1d12407274b4c4cbfd0513d44cad410c43b350a28"
  },
  "warnings": [],
  "status": "ok",
  "result": {
    "format": "tsv",
    "columns": [
      "variant_id",
      "chrom",
      "pos",
      "beta",
      "p"
    ],
    "n_records": 2,
    "text": "variant_id\tchrom\tpos\tbeta\tp\nrs1\t1\t100000\t0.12\t2e-09\nrs2\t1\t150000\t0.09\t6.8e-06\n",
    "path": null
  }
}
```
