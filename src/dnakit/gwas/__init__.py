"""Native GWAS summary-statistics and association-analysis functions.

Tables are lists of dict records, dicts of equal-length columns, delimited text/paths (.gz
allowed). Recognised column aliases (case/punctuation-insensitive): SNP/rsid -> variant_id, CHR
-> chrom, BP/POS -> pos, A1/EA/ALT -> effect_allele, A2/NEA/REF -> other_allele, BETA, SE,
P/PVAL, EAF/FREQ, N, Z, OR (converted to beta = ln OR). Every function returns the DNAKit
candidate envelope (``status``, ``result``, ``warnings``).
"""

from __future__ import annotations

from collections.abc import Mapping

from . import _impl_models as _models
from . import _impl_sumstats as _sumstats
from ._impl_common import finish as _finish

_M = "gwas"


def annotate_lead_snps(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate genome-wide-significant lead SNPs with the nearest and overlapping genes from a
    user-supplied gene table.

    input: summary-statistics table. options: genes (required: records with name/chrom/start/end
    and optional strand; no database is queried), leads (explicit variant ids; default = leads
    selected with p_threshold=5e-8 and window=500000), flank=100000 (genes listed within this
    distance), max_distance=1000000 (beyond it no nearest gene is reported). Distance is 0
    inside a gene, otherwise the distance to the nearest gene boundary.
    """

    return _finish(_M, "annotate_lead_snps", input, options, _sumstats.annotate_lead_snps)


def apply_genomic_control(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Apply genomic control: divide chi-square statistics by the inflation factor lambda_GC and
    recompute p-values.

    input: summary-statistics table with p (or z, or beta+se). lambda_GC = median(chi2) /
    0.4549364 where chi2 is the 1-df chi-square implied by each p. options: lambda_gc (override
    the estimate), only_inflated=True (never deflate: use max(1, lambda)). Adds p_gc, se_gc =
    se*sqrt(lambda) and z_gc = z/sqrt(lambda) to each record.
    """

    return _finish(_M, "apply_genomic_control", input, options, _sumstats.apply_genomic_control)


def build_credible_set(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a fine-mapping credible set from per-variant beta/se using Wakefield approximate
    Bayes factors (single causal variant).

    input: summary statistics of the variants in ONE locus (needs beta and se). log ABF =
    0.5*(ln(V/(V+W)) + W/(V+W)*z^2) with V = se^2 and W = prior_sd^2; posterior inclusion
    probability = ABF / sum(ABF) (uniform prior); the set is the smallest set of top-PIP
    variants whose PIPs sum to >= coverage. options: prior_sd=0.15, coverage=0.95, region
    ('chrom:start-end') to subset the table.
    """

    return _finish(_M, "build_credible_set", input, options, _sumstats.build_credible_set)


def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate a summary statistic of GWAS summary statistics (default: a table summary with
    genomic inflation lambda_GC).

    input: summary-statistics table. options: statistic = 'summary' (default: n, min p, lead
    variant, lambda_gc, mean/median chi2, number significant), 'lambda_gc', 'z_scores',
    'neglog10p', 'maf', 'n_significant' or 'adjusted_p' (Bonferroni p*m and Benjamini-Hochberg
    q); significance=5e-8.
    """

    return _finish(_M, "calculate", input, options, _sumstats.calculate)


def calibrate_prs_risk(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calibrate a polygenic score to absolute disease risk with a logistic recalibration model.

    input: {'score': [...], 'phenotype': [0/1, ...]} (individuals with a missing value are
    dropped, >= 10 needed). Fits risk = expit(a + b*z) with z the standardised score. options:
    prevalence (population prevalence in (0,1); shifts the intercept by logit(prevalence) -
    logit(sample prevalence) to undo case-control sampling), n_bins=10 (calibration table of
    mean predicted vs observed risk), new_scores (list of scores to convert to risk).
    """

    return _finish(_M, "calibrate_prs_risk", input, options, _models.calibrate_prs_risk)


def clump_and_threshold_prs(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Compute clumping-and-thresholding (C+T) polygenic scores from summary statistics and
    dosage genotypes.

    input: {'summary_statistics': table (variant_id, chrom, pos, beta, p), 'genotypes':
    samples x variants effect-allele dosages in [0,2] or null, 'sample_ids': [...],
    'variant_ids': [...]}. PLINK-style greedy clumping by ascending p: variants on the same
    chromosome within clump_window bp with r^2 >= clump_r2 (from the LD reference) join the
    index variant. Score = sum(beta*dosage) over index variants with p <= threshold (missing
    dosage -> variant mean). options: p_thresholds=[5e-8,1e-5,1e-3,0.05,0.5,1], clump_r2=0.1,
    clump_window=250000, ld_reference (genotype mapping; default = the target genotypes).
    """

    return _finish(_M, "clump_and_threshold_prs", input, options, _models.clump_and_threshold_prs)


def compute_locus_plot_data(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Compute the points of a regional (LocusZoom-style) association plot around a lead variant.

    input: summary-statistics table. options: lead (variant id; default = smallest p),
    flank=250000 bp each side, ld ({'variant_ids','matrix'} of correlations r) or genotypes
    (reference dosages) to colour points by r^2 with the lead (bands 0-0.2,...,0.8-1, 'lead',
    'unknown').
    """

    return _finish(_M, "compute_locus_plot_data", input, options, _sumstats.compute_locus_plot_data)


def compute_manhattan_data(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute Manhattan-plot coordinates (-log10 p and cumulative genomic x position) from
    summary statistics.

    input: summary-statistics table with chrom, pos and p. x = chromosome offset (sum of
    earlier chromosomes' maximum observed positions) + pos. options: min_neglog10p=0 (drop
    weaker points), max_points=200000. Reports the genome-wide (5e-8) and suggestive (1e-5)
    thresholds and chromosome label centres.
    """

    return _finish(_M, "compute_manhattan_data", input, options, _sumstats.compute_manhattan_data)


def compute_qq_data(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute quantile-quantile plot data (observed vs expected -log10 p) and lambda_GC.

    input: summary-statistics table with p. expected = -log10((rank-0.5)/n) for ascending-sorted
    p. options: max_points=50000 (keeps the most significant half exactly and thins the rest).
    """

    return _finish(_M, "compute_qq_data", input, options, _sumstats.compute_qq_data)


def conditional_gwas(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Run an approximate conditional association analysis from summary statistics and LD
    (COJO-style, exact OLS algebra).

    input: summary-statistics table with z (or beta+se) and ideally n. options: condition_on
    (required: 1-20 variant ids), ld ({'variant_ids','matrix'} correlation matrix r) or
    genotypes (reference dosages), n (sample size; default = median n column). For each other
    variant, the standardised effect, se, t and p from the joint model {conditioning set +
    variant} (t-test, df = n-|C|-2); marginal r is z/sqrt(n-2+z^2).
    """

    return _finish(_M, "conditional_gwas", input, options, _models.conditional_gwas)


def define_association_loci(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Define independent association loci as lead SNP +/- flank, merging overlapping windows.

    input: summary-statistics table with chrom, pos, p. options: p_threshold=5e-8, window=500000
    (lead-SNP selection distance), flank=500000 (locus half-width, clamped at position 1).
    Reports each locus' leads, strongest variant, number of variants and number of significant
    variants.
    """

    return _finish(_M, "define_association_loci", input, options, _sumstats.define_association_loci)


def estimate_genetic_correlation(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Estimate the SNP-based genetic correlation between two traits.

    input: either {'genotypes': samples x variants dosages, 'phenotype_1': [...], 'phenotype_2':
    [...]} (cross-trait Haseman-Elston regression on the VanRaden GRM) or {'trait1': table,
    'trait2': table} with z-scores and LD scores (unweighted cross-trait LD-score regression;
    alleles of trait2 are aligned to trait1). options: ld_scores ({variant_id: score}) if the
    tables have no 'ldscore' column, n1, n2, m, fix_intercept=False. rg = covariance /
    sqrt(h2_1*h2_2); undefined (None) when a heritability estimate is not positive. rg is a ratio
    estimate and is NOT bounded: it can fall outside [-1, 1] (typically with small samples); the value
    is kept unclamped and a warning is added whenever |rg| > 1. Warnings are also added when an
    individual h2 estimate is < 0 or > 1, and (summary-data path) always to state that no standard
    errors are computed.
    """

    return _finish(
        _M, "estimate_genetic_correlation", input, options, _models.estimate_genetic_correlation
    )


def estimate_snp_heritability(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Estimate SNP heritability by Haseman-Elston regression (individual data) or unweighted
    LD-score regression (summary data).

    input: {'genotypes': samples x variants dosages, 'phenotype': [...]} (method
    'haseman_elston': slope of y_i*y_j on the VanRaden GRM over pairs) or a summary-statistics
    table with z and LD scores (method 'ldsc': OLS of chi2 on LD score, h2 = slope*M/N).
    options: method='auto', ld_scores ({variant_id: score}), n, m, fix_intercept=False. No
    standard errors are produced (the summary-data path always adds a warning saying so). The estimate
    is not clamped: a warning is added when h2 < 0 or h2 > 1 (noise or misfit, not a valid heritability).
    """

    return _finish(
        _M, "estimate_snp_heritability", input, options, _models.estimate_snp_heritability
    )


def evaluate_prs(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Evaluate a polygenic score against a phenotype: R^2 / slope for quantitative traits, AUC /
    odds ratio for binary traits.

    input: {'score': [...], 'phenotype': [...]}. options: trait='auto' | 'quantitative' |
    'binary'. Quantitative: Pearson r, r^2, slope per SD of score with t-test p. Binary: AUC
    (Mann-Whitney), odds ratio per SD from logistic regression with Wald 95% CI and p,
    Nagelkerke R^2.
    """

    return _finish(_M, "evaluate_prs", input, options, _models.evaluate_prs)


def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Serialise GWAS summary statistics to another file format (tsv, csv, json, jsonl, bed,
    ldsc, metal).

    input: summary-statistics table. options: format='tsv', columns (subset for tsv/csv),
    path (write the text there; .gz is gzip-compressed; existing files are not overwritten
    unless overwrite=True). 'bed' holds -log10 p as score, 'ldsc' the SNP/A1/A2/N/Z munged
    layout, 'metal' the METAL input layout.
    """

    return _finish(_M, "export", input, options, _sumstats.export_table)


def filter(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter GWAS summary-statistics rows by p-value, MAF, sample size, INFO, chromosome, region
    or variant type.

    input: summary-statistics table. options: max_p, min_maf, min_n, min_info, chromosomes
    (list), region ('chrom:start-end'), snp_only=False, remove_palindromic=False,
    require_p=False. Returns kept records and the number removed per reason.
    """

    return _finish(_M, "filter", input, options, _sumstats.filter_rows)


def fixed_effect_meta_analysis(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Combine per-variant effects across studies with inverse-variance-weighted fixed-effect
    meta-analysis.

    input: a list or dict of study tables (each with variant_id or chrom+pos, beta, se, alleles)
    or one table with a 'study' column. Alleles are aligned to the first study (effect flipped
    when swapped; mismatches skipped). Returns beta, se, z, p, Cochran's Q (with p) and I^2.
    options: min_studies=2, study_column='study'.
    """

    return _finish(
        _M, "fixed_effect_meta_analysis", input, options, _sumstats.fixed_effect_meta_analysis
    )


def harmonize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Harmonise effect alleles of a summary-statistics table to a reference table (swap /
    strand-flip / palindrome handling).

    input: {'target': table, 'reference': table} or a table with options['reference']. Matches
    variants by position or id (match_on='auto'), flips beta, z and eaf when the effect allele
    is the reference's other allele, complements alleles on the opposite strand, and drops
    allele mismatches and ambiguous palindromic (A/T, C/G) SNPs whose frequency is within
    palindrome_eaf_tolerance=0.08 of 0.5 (palindromes='drop_ambiguous', or
    'assume_forward_strand' to keep them).
    """

    return _finish(_M, "harmonize", input, options, _sumstats.harmonize)


def joint_gwas_analysis(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Run a joint multi-SNP association model from summary statistics and LD (COJO-joint style).

    input: summary-statistics table with z (or beta+se) and ideally n. options: ld
    ({'variant_ids','matrix'} correlation matrix r) or genotypes (reference dosages), snps
    (variant ids to model jointly, max 200; default all in the LD matrix), n. Computes the exact
    multiple-regression coefficients of the standardised phenotype from marginal correlations:
    beta = R^-1 r, se^2 = (1-R^2)/(n-k-1) * diag(R^-1), t-test p.
    """

    return _finish(_M, "joint_gwas_analysis", input, options, _models.joint_gwas_analysis)


def merge_nearby_significant_hits(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Merge nearby genome-wide-significant variants into hit regions by single-linkage on
    base-pair distance.

    input: summary-statistics table with chrom, pos, p. options: p_threshold=5e-8,
    distance=500000. Consecutive significant variants on one chromosome at most `distance` bp
    apart share a region (start, end, lead variant, count).
    """

    return _finish(
        _M, "merge_nearby_significant_hits", input, options, _sumstats.merge_nearby_significant_hits
    )


def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw a Manhattan or QQ plot of GWAS p-values as an SVG string.

    input: summary-statistics table. options: kind='manhattan' (needs chrom, pos, p) or 'qq',
    width=900, height=450, max_points=20000, path (write the SVG there; not overwritten unless
    overwrite=True).
    """

    return _finish(_M, "plot", input, options, _sumstats.plot)


def random_effect_meta_analysis(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Combine per-variant effects across studies with DerSimonian-Laird random-effects
    meta-analysis.

    input and options as fixed_effect_meta_analysis. tau^2 = max(0, (Q-df)/(S1 - S2/S1)) with S1
    = sum(w), S2 = sum(w^2), w = 1/se^2; pooled effect uses weights 1/(se^2+tau^2) with a
    normal-theory p-value.
    """

    return _finish(
        _M, "random_effect_meta_analysis", input, options, _sumstats.random_effect_meta_analysis
    )


def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Read a GWAS summary-statistics table (path, .gz file, delimited text, records or columns)
    into records.

    Delimiter is detected (tab, comma, semicolon, whitespace); NA/empty cells become null and
    numeric text becomes numbers. options: delimiter, max_rows=500000, standardize=False (rename
    to the canonical columns).
    """

    return _finish(_M, "read", input, options, _sumstats.read_table)


def run(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Run a native single-variant association scan (linear or logistic regression) on dosage
    genotypes.

    input: {'genotypes': samples x variants dosages in [0,2] or null, 'phenotype': [...],
    'covariates': samples x k (optional), 'variant_ids', 'chrom', 'pos', 'alt', 'ref',
    'sample_ids'}. Per variant: OLS t-test (quantitative trait, df = n-k-2) or logistic Wald
    z-test (0/1 trait) on the dosage coefficient, with intercept and covariates; samples with
    missing genotype or phenotype are dropped per variant. options: model='auto' | 'linear' |
    'logistic', maf_min=0, backend='native'. External programs (PLINK, REGENIE, GEMMA, GCTA) are
    not driven: a requested external backend reports not_implemented.
    """

    return _finish(_M, "run", input, options, _models.run_association)


def select_lead_snps(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Select independent lead SNPs by greedy distance-based pruning of genome-wide-significant
    variants.

    input: summary-statistics table with chrom, pos, p. Variants with p <= p_threshold are taken
    in order of ascending p; each becomes a lead and removes other candidates within `window` bp
    on the same chromosome. options: p_threshold=5e-8, window=500000, max_leads.
    """

    return _finish(_M, "select_lead_snps", input, options, _sumstats.select_lead_snps)


def simulate_gwas_data(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Simulate genotypes, a quantitative phenotype and the resulting association summary
    statistics (deterministic by seed).

    input: None or a mapping of the options below (overrides options). Dosage ~ Binomial(2, f_j)
    with f_j uniform in [maf_min, maf_max]; n_causal variants get effects ~ N(0,
    heritability/n_causal) on standardised dosage; phenotype = genetic value + N(0,
    1-heritability). options: n_samples=500, n_variants=200, n_causal=5, heritability=0.5,
    seed=0, maf_min=0.05, maf_max=0.5, n_chromosomes=22. Returns genotypes, phenotype, true
    effects and the linear-regression scan.
    """

    return _finish(_M, "simulate_gwas_data", input, options, _models.simulate_gwas_data)


def standardize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Standardise GWAS summary-statistics columns to the canonical names and derive missing
    statistics.

    Canonical columns: variant_id, chrom (no 'chr' prefix; 23->X, 24->Y, 26->MT), pos,
    effect_allele, other_allele (upper case), beta, se, p, eaf, n, z. Derives beta = ln(OR), z =
    beta/se, p from z (two-sided normal), se from beta and p, and variant_id =
    chrom:pos:other:effect when absent. options: max_rows=500000, delimiter.
    """

    return _finish(_M, "standardize", input, options, _sumstats.standardize)


def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate a GWAS summary-statistics table: required columns, ranges, alleles and duplicate
    variants.

    Reports errors (no id/position, no statistic, p outside [0,1], se <= 0, eaf outside [0,1], n
    <= 0, non-finite values, identical alleles) and warnings (non-ACGT alleles, duplicate ids,
    missing allele/eaf columns). Returns valid=False rather than raising when the content is
    wrong; unreadable tables raise.
    """

    return _finish(_M, "validate", input, options, _sumstats.validate)


def write(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Write GWAS summary statistics as a delimited table (tsv/csv) with canonical column names.

    options: path (output file; .gz is gzip-compressed; the text is always returned; existing
    files are not overwritten unless overwrite=True), format='tsv' ('csv' for .csv paths),
    columns (subset).
    """

    return _finish(_M, "write", input, options, _sumstats.write_table)


__all__ = [
    "annotate_lead_snps",
    "apply_genomic_control",
    "build_credible_set",
    "calculate",
    "calibrate_prs_risk",
    "clump_and_threshold_prs",
    "compute_locus_plot_data",
    "compute_manhattan_data",
    "compute_qq_data",
    "conditional_gwas",
    "define_association_loci",
    "estimate_genetic_correlation",
    "estimate_snp_heritability",
    "evaluate_prs",
    "export",
    "filter",
    "fixed_effect_meta_analysis",
    "harmonize",
    "joint_gwas_analysis",
    "merge_nearby_significant_hits",
    "plot",
    "random_effect_meta_analysis",
    "read",
    "run",
    "select_lead_snps",
    "simulate_gwas_data",
    "standardize",
    "validate",
    "write",
]

__candidate_adapter_operations__ = set(
    [
        "annotate_lead_snps",
        "apply_genomic_control",
        "build_credible_set",
        "calculate",
        "calibrate_prs_risk",
        "clump_and_threshold_prs",
        "compute_locus_plot_data",
        "compute_manhattan_data",
        "compute_qq_data",
        "conditional_gwas",
        "define_association_loci",
        "estimate_genetic_correlation",
        "estimate_snp_heritability",
        "evaluate_prs",
        "export",
        "filter",
        "fixed_effect_meta_analysis",
        "harmonize",
        "joint_gwas_analysis",
        "merge_nearby_significant_hits",
        "plot",
        "random_effect_meta_analysis",
        "read",
        "run",
        "select_lead_snps",
        "simulate_gwas_data",
        "standardize",
        "validate",
        "write",
    ]
)


from dnakit._verbosity import register_public as _register_public_verbose

_register_public_verbose(globals())
