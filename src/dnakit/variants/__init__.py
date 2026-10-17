"""Local candidate adapters for the variants capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok, run_candidate_operation
from dnakit._verbosity import register_public as _register_public_verbose

from .external import (
    annotate_vcf_ensembl_vep,
    annotate_vcf_vcfanno,
    call_cutesv,
    call_deepvariant,
    call_delly,
    call_freebayes,
    call_gatk_haplotypecaller,
    call_sniffles2,
    filter_vcf_vcftools,
    load_gemini_database,
    merge_structural_variants_survivor,
    query_gemini_database,
)
from .statistics import (
    calculate_allele_frequency,
    calculate_heterozygosity,
    calculate_missingness,
    calculate_transition_transversion,
    calculate_variant_density,
    calculate_minor_allele_frequency as _calculate_minor_allele_frequency_local,
    compare as _compare_local,
    count as _count_local,
    filter_variants_by_allele_frequency as _filter_variants_by_allele_frequency_local,
    filter_variants_by_depth as _filter_variants_by_depth_local,
    filter_variants_by_genotype_quality as _filter_variants_by_genotype_quality_local,
    filter_variants_by_missingness as _filter_variants_by_missingness_local,
    filter_variants_by_quality as _filter_variants_by_quality_local,
    parse_vcf_records as _parse_vcf_records_local,
    split_multiallelic as _split_multiallelic_local,
)
from .backends import (
    calculate_fst,
    calculate_nucleotide_diversity,
    calculate_observed_heterozygosity,
    calculate_pca,
    calculate_site_frequency_spectrum,
    calculate_tajimas_d,
    calculate_watterson_theta,
    calculate_windowed_ld,
    fetch_vcf_region,
    read_vcf,
)

from . import _impl_analysis as _analysis
from . import _impl_prioritize as _prior
from . import _impl_vcf as _vcf


def _done(operation: str, input: object, options: Mapping[str, object], result: dict[str, object]) -> dict[str, object]:
    """Wrap a computed result in the candidate envelope (``warnings`` move to the envelope)."""

    warnings = list(result.pop("warnings", []) or [])
    return candidate_ok("variants", operation, input, result, options=options, warnings=warnings)


def _vcf_result(vs: "_vcf.VariantSet", options: Mapping[str, object], **extra: object) -> dict[str, object]:
    """Result dict for operations that return a transformed call set; writes ``options['output']`` if given."""

    text = _vcf.to_vcf_text(vs)
    result: dict[str, object] = {"n_records": len(vs.records), "samples": list(vs.samples), "vcf": text,
                                 "output_path": _vcf.write_output(options, text)}
    if _vcf.option_bool(options, "include_records", True):
        result["records"] = [_vcf.record_view(rec, vs.samples) for rec in vs.records]
    result.update(extra)
    return result


def _add_hgvs(vs: "_vcf.VariantSet", options: Mapping[str, object], reference: "_vcf.Reference | None") -> None:
    accession = options.get("hgvs_accession")
    if accession is None:
        return
    if not isinstance(accession, (str, Mapping)):
        raise _vcf.ConfigurationError("options.hgvs_accession must be a string or {contig: accession}.",
                                      code="VARIANTS_BAD_OPTION", hint="Example: 'NC_000001.11'.")
    for rec in vs.records:
        acc = accession.get(rec["chrom"]) if isinstance(accession, Mapping) else accession
        if acc:
            rec["hgvs_g"] = [_vcf.build_hgvs_genomic(acc, rec["pos"], rec["ref"], alt, reference, rec["chrom"])
                             for alt in rec["alts"] if not _vcf.is_symbolic(alt)]


def _normalize_impl(input: object, options: Mapping[str, object], *, require_reference: bool) -> dict[str, object]:
    vs = _vcf.load_variant_set(input)
    reference = _vcf.reference_from_options(options)
    if require_reference and reference is None:
        raise _vcf.ConfigurationError("left-alignment needs a reference sequence.", code="VARIANTS_REFERENCE_REQUIRED",
                                      hint="Pass options.reference (a sequence, FASTA text/path or {contig: sequence}).")
    out, report = _vcf.normalize_set(
        vs, reference, check_ref=_vcf.option_choice(options, "check_ref", "warn", ("warn", "error")),
        split=_vcf.option_bool(options, "split_multiallelic", False), sort=_vcf.option_bool(options, "sort", True))
    _add_hgvs(out, options, reference)
    warnings = report.pop("warnings")
    return _vcf_result(out, options, warnings=warnings, **report)


def _compare_vcfs_impl(input: object, options: Mapping[str, object]) -> dict[str, object]:
    if isinstance(input, Mapping) and "left" in input and "right" in input:
        left_in, right_in = input["left"], input["right"]
    elif options.get("other") is not None:
        left_in, right_in = input, options["other"]
    else:
        raise _vcf.InputFormatError("comparison needs two call sets.", code="VARIANTS_BAD_INPUT",
                                    hint="Pass {'left': ..., 'right': ...} or options.other.")
    left, right = _vcf.load_variant_set(left_in, what="left"), _vcf.load_variant_set(right_in, what="right")
    reference = _vcf.reference_from_options(options)
    warnings: list[str] = []
    if reference is not None:
        check = _vcf.option_choice(options, "check_ref", "warn", ("warn", "error"))
        left, rep_l = _vcf.normalize_set(left, reference, check_ref=check)
        right, rep_r = _vcf.normalize_set(right, reference, check_ref=check)
        warnings = [*rep_l["warnings"], *rep_r["warnings"]]
    result = _vcf.compare_sets(left, right, match=_vcf.option_choice(options, "match", "exact", ("exact", "position")))
    result["normalized_with_reference"] = reference is not None
    result["warnings"] = warnings
    return result


def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate variants with functional, population and clinical information (not implemented).

    Needs transcript models, population frequencies and clinical databases; use the specific annotate_* functions or annotate_vcf_ensembl_vep / annotate_vcf_vcfanno.

    Currently returns status 'not_implemented': annotation needs transcript, population-frequency and clinical reference data that DNAKit does not bundle.
    """

    return candidate_not_implemented('variants', 'annotate', input, options=options,
        requires='annotation databases (Ensembl/VEP cache, gnomAD, ClinVar) or the vep / vcfanno programs',
        message='annotation needs transcript, population-frequency and clinical reference data that DNAKit does not bundle.')

def annotate_coding_effect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Predict the protein-level effect (synonymous, missense, nonsense, start/stop loss, frameshift, in-frame indel) of variants on a user-supplied coding sequence.

    input: {'cds': 'ATG...', 'variants': [{'cds_position': 5, 'ref': 'C', 'alt': 'T'}]} or 'c.5C>T' strings; options.cds may supply the CDS instead.
    options: cds (coding sequence on the coding strand, ACGT).
    Positions are 1-based CDS coordinates, alleles are on the coding strand, the standard genetic code is used and HGVS protein strings are not 3'-normalised. No transcript model or genome is used; for genomic variants use annotate_vcf_ensembl_vep.
    """

    options = dict(options or {})
    return _done('annotate_coding_effect', input, options, _analysis.coding_effect(input, options))

def annotate_conservation_score(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate variants with phyloP/phastCons/GERP conservation scores (not implemented).

    Needs a genome-wide conservation track (bigWig/bedGraph).

    Currently returns status 'not_implemented': conservation scoring needs a genome-wide conservation track that is not available.
    """

    return candidate_not_implemented('variants', 'annotate_conservation_score', input, options=options,
        requires='per-base conservation track (phyloP, phastCons or GERP)',
        message='conservation scoring needs a genome-wide conservation track that is not available.')

def annotate_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate mitochondrial variants (not implemented).

    Needs MITOMAP/gnomAD-mt/HmtVar style mitochondrial annotation resources.

    Currently returns status 'not_implemented': mitochondrial annotation needs reference annotation databases that are not available.
    """

    return candidate_not_implemented('variants', 'annotate_mt', input, options=options,
        requires='mitochondrial annotation resources (MITOMAP, gnomAD-mt, HmtVar)',
        message='mitochondrial annotation needs reference annotation databases that are not available.')

def annotate_pathogenicity_score(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate variants with a pathogenicity prediction score (not implemented).

    Needs a trained predictor (CADD, REVEL, AlphaMissense, ...) or its precomputed score files.

    Currently returns status 'not_implemented': pathogenicity scores come from trained models or precomputed score files that are not available.
    """

    return candidate_not_implemented('variants', 'annotate_pathogenicity_score', input, options=options,
        requires='trained pathogenicity predictor or precomputed score tables (CADD/REVEL/AlphaMissense)',
        message='pathogenicity scores come from trained models or precomputed score files that are not available.')

def annotate_population_frequency(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate variants with population allele frequencies (not implemented).

    Needs gnomAD/1000 Genomes/dbSNP frequency data (a local file or a network query).

    Currently returns status 'not_implemented': population frequencies need a frequency database or network access, neither of which is used here.
    """

    return candidate_not_implemented('variants', 'annotate_population_frequency', input, options=options,
        requires='population frequency database (gnomAD, 1000 Genomes, dbSNP)',
        message='population frequencies need a frequency database or network access, neither of which is used here.')

def annotate_regulatory_effect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate variants with regulatory-element effects (not implemented).

    Needs regulatory build annotations (Ensembl Regulatory Build, ENCODE cCREs) or a sequence-to-activity model.

    Currently returns status 'not_implemented': regulatory annotation needs regulatory-region tracks or trained models that are not available.
    """

    return candidate_not_implemented('variants', 'annotate_regulatory_effect', input, options=options,
        requires='regulatory annotation track or trained regulatory model',
        message='regulatory annotation needs regulatory-region tracks or trained models that are not available.')

def annotate_splice_effect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate variants with predicted splicing effects (not implemented).

    Needs exon/transcript models plus a trained splice predictor (SpliceAI, MaxEntScan data).

    Currently returns status 'not_implemented': splice-effect prediction needs transcript models and a trained predictor that are not available.
    """

    return candidate_not_implemented('variants', 'annotate_splice_effect', input, options=options,
        requires='transcript annotation and a splice-effect predictor (SpliceAI)',
        message='splice-effect prediction needs transcript models and a trained predictor that are not available.')

def annotate_variant_consequence(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate Sequence Ontology consequences of variants on transcripts (not implemented).

    Needs transcript models (GTF/GFF + genome FASTA). For a single supplied CDS use annotate_coding_effect.

    Currently returns status 'not_implemented': consequence calling needs transcript models and a reference genome; only annotate_coding_effect works from a user-supplied CDS.
    """

    return candidate_not_implemented('variants', 'annotate_variant_consequence', input, options=options,
        requires='transcript annotation (GTF/GFF) and genome FASTA, e.g. via Ensembl VEP',
        message='consequence calling needs transcript models and a reference genome; only annotate_coding_effect works from a user-supplied CDS.')

def calculate_copy_number_profile(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Convert per-bin sequencing depth into a copy-number profile (depth ratio to the median, copy number, merged segments).

    input: a list of bins [{'chrom','start','end','depth'[,'control_depth']}] (or [chrom,start,end,depth]) or {'bins': [...]}.
    options: purity (0.01-1, default 1.0), ploidy (default 2.0), median_exclude_chroms (default X, Y, M, MT).
    Ratio = depth / median depth (divided by the control's ratio when control_depth is given); copy number = (ratio*(purity*ploidy+2*(1-purity)) - 2*(1-purity))/purity.
    Adjacent bins with the same rounded copy number are merged into segments; no statistical segmentation (CBS/HMM) is performed.
    """

    options = dict(options or {})
    return _done('calculate_copy_number_profile', input, options, _analysis.copy_number_profile(input, options))

def calculate_genotype_likelihood(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute genotype likelihoods (log10), Phred-scaled PL values and GQ for a biallelic site from read evidence under a binomial error model.

    input: {'ref_count': int, 'alt_count': int} or {'reads': [{'allele': 'ref'|'alt', 'base_quality': q}]} (or a list of those), or a VCF/records with FORMAT/AD per sample.
    options: error_rate (default 0.01, used when base qualities are absent), ploidy (default 2).
    Each read supports ALT with probability (k/m)(1-e)+(1-k/m)e under a genotype with k ALT copies of m; PL = -10*(log10 L - max log10 L); GQ = second-lowest PL capped at 99.
    """

    options = dict(options or {})
    return _done('calculate_genotype_likelihood', input, options, _analysis.genotype_likelihood(input, options))

def calculate_haplotype_frequency(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate haplotype frequencies at two or more loci: direct counting of phased haplotypes or EM estimation from unphased genotypes.

    input: {'haplotypes': ['010', '011', ...]} (phased), {'genotypes': [['0/1','1/1'], ...]} or a list of per-individual GT lists, or a VCF/records with sample genotypes (loci = biallelic records).
    options: use_phase (default True, VCF input: count haplotypes directly when every genotype is phased), max_iterations (1000), tolerance (1e-10).
    EM follows Excoffier & Slatkin (1995) under Hardy-Weinberg with a uniform start; individuals with missing genotypes are excluded (up to 12 loci, 10 heterozygous loci per individual).
    """

    options = dict(options or {})
    return _done('calculate_haplotype_frequency', input, options, _analysis.haplotype_frequency(input, options))

def calculate_hwe_pvalue(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Test Hardy-Weinberg equilibrium with the exact test of Wigginton, Cutler & Abecasis (2005).

    input: genotype counts [hom_ref, het, hom_alt] (or {'counts': [...]}, {'hom_ref','het','hom_alt'}), or a VCF/records with sample genotypes (one p-value per biallelic diploid site).
    options: mid_p (bool, default False; subtracts half the probability of the observed table).
    The two-sided p-value sums the probabilities of all heterozygote counts (given the allele counts) that are not more likely than the observed one.
    """

    options = dict(options or {})
    return _done('calculate_hwe_pvalue', input, options, _analysis.hwe_pvalue(input, options))

def calculate_kinship_king(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate pairwise kinship coefficients between samples with the KING-robust estimator (a genotype-based estimate of the kinship coefficient, not IBD segments).

    input: a multi-sample VCF (text/path) or records with diploid genotypes (at most 150 samples; biallelic sites are used).
    options: none.
    kinship = (N_AaAa - 2*N_AAaa) / (N_Aa(i) + N_Aa(j)) over sites called in both samples (Manichaikul et al. 2010); pi_hat_approx = 2*kinship; also IBS0 fraction and a degree label by kinship cut-offs (duplicate/monozygotic twin > 0.3535, 1st degree > 0.1768, 2nd degree > 0.0884, 3rd degree > 0.0442, else unrelated; the KING boundaries 2^-1.5, 2^-2.5, 2^-3.5 and 2^-4.5 given to four decimals; a value equal to a boundary falls in the lower class). It does not detect or report IBD segments; pi_hat_approx is only 2*kinship, not a PLINK PI_HAT.
    """

    options = dict(options or {})
    return _done('calculate_kinship_king', input, options, _analysis.identity_by_descent(input, options))

def calculate_ld_matrix(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compute a pairwise linkage-disequilibrium matrix (r2, r or D') between biallelic variants.

    input: VCF text/path or records with sample genotypes (at least two samples and two biallelic variants).
    options: statistic ('r2' default, 'r', 'd_prime'), phase ('auto' default, 'phased', 'unphased'), max_variants (default 200, max 1000).
    Phased genotypes give haplotype-based D, D' and r2; unphased genotypes give the Pearson correlation of ALT dosages (Rogers-Huff r), D' is unavailable.
    """

    options = dict(options or {})
    return _done('calculate_ld_matrix', input, options, _analysis.ld_matrix(input, options))

def calculate_runs_of_homozygosity(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Find runs of homozygosity in one sample from consecutive homozygous genotypes.

    input: VCF text/path or records with genotypes.
    options: sample (default first), min_snps (25), min_length_bp (1000000), max_heterozygotes (1), max_missing (2), max_gap_bp (1000000).
    Greedy left-to-right scan per chromosome: a run starts at a homozygous call and extends while heterozygous/missing calls stay within the allowed counts and gaps stay below max_gap_bp; trailing heterozygous/missing calls are trimmed. Runs shorter than min_snps or min_length_bp are discarded. This is a consecutive-SNP heuristic, not an HMM.
    """

    options = dict(options or {})
    return _done('calculate_runs_of_homozygosity', input, options, _analysis.runs_of_homozygosity(input, options))

def call(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Call variants from sequencing reads (not implemented).

    Needs a variant caller; use call_freebayes, call_gatk_haplotypecaller, call_deepvariant, call_delly, call_cutesv or call_sniffles2, which run the external programs when they are installed.

    Currently returns status 'not_implemented': variant calling needs an external caller and aligned reads; use the call_<caller> functions.
    """

    return candidate_not_implemented('variants', 'call', input, options=options,
        requires='variant caller binary (bcftools/FreeBayes/GATK/DeepVariant) and an aligned BAM/CRAM with a reference',
        message='variant calling needs an external caller and aligned reads; use the call_<caller> functions.')

def call_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Call mitochondrial variants from reads (not implemented).

    Needs a mitochondrial-aware caller (GATK Mutect2 mitochondria mode, mutserve) and aligned reads.

    Currently returns status 'not_implemented': mitochondrial variant calling needs an external caller and aligned reads.
    """

    return candidate_not_implemented('variants', 'call_mt', input, options=options,
        requires='mitochondrial variant caller (Mutect2 --mitochondria-mode, mutserve) and aligned reads',
        message='mitochondrial variant calling needs an external caller and aligned reads.')

def compare_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare two variant-by-sample (matrix-table style) call sets: shared and private alleles, precision, recall, F1 and genotype concordance.

    input: {'left': <VCF|records>, 'right': <VCF|records>} (left is the truth/reference set), or one call set with the other in options.other.
    options: other, match ('exact' default: CHROM/POS/REF/ALT, or 'position'), reference (+ reference_start) to normalise both sets first, check_ref.
    Genotype concordance is computed per sample name present in both sets over shared sites (alleles compared by sequence, order-insensitive).
    """

    options = dict(options or {})
    return _done('compare_mt', input, options, _compare_vcfs_impl(input, options))

def compare_vcfs(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare two VCFs / variant sets: shared and private alleles, precision, recall, F1 and genotype concordance.

    input: {'left': <VCF|records>, 'right': <VCF|records>} (left is the truth/reference set), or one call set with the other in options.other.
    options: other, match ('exact' default: CHROM/POS/REF/ALT, or 'position'), reference (+ reference_start) to normalise both sets first, check_ref.
    Genotype concordance is computed per sample name present in both sets over shared sites (alleles compared by sequence, order-insensitive).
    """

    options = dict(options or {})
    return _done('compare_vcfs', input, options, _compare_vcfs_impl(input, options))

def construct_haplotypes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Construct per-sample haplotypes (allele lists and, with a reference, sequences) from phased genotypes.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples). Genotypes should be phased (0|1).
    options: sample/samples, on_unphased ('error' default or 'skip'), reference (+ reference_start) and window_start/window_end to also build haplotype sequences.
    The first GT allele is haplotype 1. Unphased heterozygotes are an error (or skipped); this does not phase genotypes (see phase_genotypes).
    """

    options = dict(options or {})
    return _done('construct_haplotypes', input, options, _analysis.construct_haplotypes(input, options))

def count_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation('variants', 'count_mt', input, options=options)

def detect_breakpoints(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Detect structural-variant breakpoints from alignments (not implemented).

    Needs split/discordant read evidence from a BAM/CRAM and a breakpoint caller (DELLY, Manta, GRIDSS).

    Currently returns status 'not_implemented': breakpoint detection needs read alignments and an SV caller; call_delly, call_cutesv and call_sniffles2 run such tools.
    """

    return candidate_not_implemented('variants', 'detect_breakpoints', input, options=options,
        requires='SV caller (DELLY/Manta/GRIDSS) and aligned reads',
        message='breakpoint detection needs read alignments and an SV caller; call_delly, call_cutesv and call_sniffles2 run such tools.')

def detect_low_frequency_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Detect low-frequency variants from per-sample allelic depths with a binomial error-rate test.

    input: VCF text/path or records with FORMAT/AD (ref,alt) per sample.
    options: min_vaf (0.01), max_vaf (0.2), min_depth (100), min_alt_count (3), error_rate (0.001), alpha (0.01), samples/sample.
    A call is reported when min_vaf <= alt/depth <= max_vaf, depth >= min_depth, alt >= min_alt_count and P(X >= alt | depth, error_rate) <= alpha. Heuristic: no strand, context or UMI error model.
    """

    options = dict(options or {})
    return _done('detect_low_frequency_variants', input, options, _analysis.low_frequency_variants(input, options))

def detect_mobile_element_insertions(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Select mobile-element insertion records (Alu, LINE-1, SVA, HERV) from a structural-variant VCF.

    input: SV VCF text/path or records with symbolic ALT such as <INS:ME:ALU> or INFO SVTYPE/MEINFO.
    options: families (optional list such as ['Alu', 'LINE-1']).
    This identifies MEI calls already made by a caller; it does not search alignments for insertions.
    """

    options = dict(options or {})
    return _done('detect_mobile_element_insertions', input, options, _analysis.mobile_element_insertions(input, options))

def detect_mosaic_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Detect candidate mosaic variants: ALT fractions significantly above sequencing error and significantly below a germline 0.5.

    input: VCF text/path or records with FORMAT/AD (ref,alt) per sample.
    options: min_vaf (0.03), max_vaf (0.4), min_depth (20), min_alt_count (3), error_rate (0.002), alpha (0.01), samples/sample.
    Heuristic: needs both the one-sided binomial test against error (upper tail) and against p=0.5 (lower tail) at level alpha; no population or strand-bias filters.
    """

    options = dict(options or {})
    return _done('detect_mosaic_variants', input, options, _analysis.mosaic_variants(input, options))

def detect_repeat_expansions(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Find the longest uninterrupted tandem repeat of a given unit in DNA sequences and flag expansions above a threshold.

    input: a DNA string, FASTA text, a list of strings or {name: sequence}.
    options: motif (required, 1-12 bases, or a list), expansion_threshold (units; optional), both_strands (default False).
    Exact repeats only (interruptions end a run); this reads a locus sequence, it does not genotype repeats from reads like ExpansionHunter.
    """

    options = dict(options or {})
    return _done('detect_repeat_expansions', input, options, _analysis.repeat_expansions(input, options))

def filter(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter variant records (site and genotype level) with quality, depth, missingness, allele-frequency, HWE, type and region criteria.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options (all optional): pass_only, min_qual, min_depth, max_depth (site DP), min_genotype_quality and min_sample_depth (calls below are set to ./.), max_missing (fraction), min_af, max_af (ALT frequency from genotypes, else INFO/AF), min_hwe_pvalue, variant_types (SNV/MNV/indel/SV), chroms, regions ('chr:start-end'), mode ('remove' default, or 'flag' to write the failing reasons into FILTER), output (VCF path), include_records, overwrite.
    Missing QUAL/depth fails the corresponding threshold. Result: filtered records, VCF text and per-reason counts.
    """

    options = dict(options or {})
    vs = _vcf.load_variant_set(input)
    out, report = _vcf.filter_set(vs, options)
    return _done('filter', input, options, _vcf_result(out, options, **report))

def filter_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Filter a variant-by-sample (matrix-table style) call set with quality, depth, missingness, allele-frequency, HWE, type and region criteria.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options (all optional): pass_only, min_qual, min_depth, max_depth (site DP), min_genotype_quality and min_sample_depth (calls below are set to ./.), max_missing (fraction), min_af, max_af (ALT frequency from genotypes, else INFO/AF), min_hwe_pvalue, variant_types (SNV/MNV/indel/SV), chroms, regions ('chr:start-end'), mode ('remove' default, or 'flag' to write the failing reasons into FILTER), output (VCF path), include_records, overwrite.
    Missing QUAL/depth fails the corresponding threshold. Result: filtered records, VCF text and per-reason counts.
    """

    options = dict(options or {})
    vs = _vcf.load_variant_set(input)
    out, report = _vcf.filter_set(vs, options)
    return _done('filter_mt', input, options, _vcf_result(out, options, **report))

def genotype_structural_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Genotype known structural variants in a sample from its reads (not implemented).

    Needs an SV genotyper (SVTyper, Paragraph, GraphTyper) and aligned reads.

    Currently returns status 'not_implemented': SV genotyping needs read alignments and an external genotyper.
    """

    return candidate_not_implemented('variants', 'genotype_structural_variants', input, options=options,
        requires='SV genotyper (SVTyper/Paragraph) and aligned reads',
        message='SV genotyping needs read alignments and an external genotyper.')

def fill_missing_genotypes_with_mode(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Fill missing genotype calls with the most frequent called genotype at the same site (naive modal imputation).

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: output (optional path to write the VCF), overwrite (default False).
    Naive and not LD-aware: it is not equivalent to reference-panel imputation (Beagle, IMPUTE5); imputed calls get FORMAT/IMPUTED=1.
    """

    options = dict(options or {})
    return _done('fill_missing_genotypes_with_mode', input, options, _analysis.impute_missing(input, options))

def left_align(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Left-align and trim indels in variant records against a reference sequence.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: reference (plain sequence, FASTA text/path or {contig: sequence}), reference_start (1-based position of a plain sequence's first base, default 1), check_ref ('warn' default: records whose REF does not match are left unchanged; 'error'), split_multiallelic (False), sort (True), hgvs_accession (accession string or {contig: accession}) to add an HGVS g. string per ALT, output (VCF path), include_records, overwrite.
    Algorithm of Tan, Abecasis & Kang (2015): drop the shared last base, extend left by the reference base when an allele becomes empty, then drop shared leading bases. Multi-allelic records are handled jointly. Required: options.reference.
    """

    options = dict(options or {})
    return _done('left_align', input, options, _normalize_impl(input, options, require_reference=True))

def left_align_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Left-align and trim indels in a variant-by-sample (matrix-table style) call set against a reference sequence.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: reference (plain sequence, FASTA text/path or {contig: sequence}), reference_start (1-based position of a plain sequence's first base, default 1), check_ref ('warn' default: records whose REF does not match are left unchanged; 'error'), split_multiallelic (False), sort (True), hgvs_accession (accession string or {contig: accession}) to add an HGVS g. string per ALT, output (VCF path), include_records, overwrite.
    Algorithm of Tan, Abecasis & Kang (2015): drop the shared last base, extend left by the reference base when an allele becomes empty, then drop shared leading bases. Multi-allelic records are handled jointly. Required: options.reference.
    """

    options = dict(options or {})
    return _done('left_align_mt', input, options, _normalize_impl(input, options, require_reference=True))

def left_align_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Left-align and trim indels in variant records against a reference sequence.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: reference (plain sequence, FASTA text/path or {contig: sequence}), reference_start (1-based position of a plain sequence's first base, default 1), check_ref ('warn' default: records whose REF does not match are left unchanged; 'error'), split_multiallelic (False), sort (True), hgvs_accession (accession string or {contig: accession}) to add an HGVS g. string per ALT, output (VCF path), include_records, overwrite.
    Algorithm of Tan, Abecasis & Kang (2015): drop the shared last base, extend left by the reference base when an allele becomes empty, then drop shared leading bases. Multi-allelic records are handled jointly. Required: options.reference.
    """

    options = dict(options or {})
    return _done('left_align_variants', input, options, _normalize_impl(input, options, require_reference=True))

def merge(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge several VCFs / call sets, combining samples (default) or concatenating records.

    input: a list of at least two call sets (VCF text, VCF paths or {'records': [...]}) or {'inputs': [...]}.
    options: mode ('samples' default, like bcftools merge -m none: identical CHROM/POS/REF/ALT sites are combined, sample names must be unique, absent samples get ./.; 'records', like bcftools concat: identical sample lists, duplicate sites kept once), output (VCF path), include_records, overwrite.
    Merged QUAL = max, FILTER = PASS only if every contributing record passes (else the union of failing labels), INFO = first record carrying the key. Output is sorted.
    """

    options = dict(options or {})
    mode = _vcf.option_choice(options, 'mode', 'samples', ('samples', 'records'))
    merged = _vcf.merge_sets(_vcf.load_many(input), mode)
    return _done('merge', input, options, _vcf_result(merged, options, mode=mode))

def merge_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge several variant-by-sample (matrix-table style) call sets, combining samples or concatenating records.

    input: a list of at least two call sets (VCF text, VCF paths or {'records': [...]}) or {'inputs': [...]}.
    options: mode ('samples' default, like bcftools merge -m none: identical CHROM/POS/REF/ALT sites are combined, sample names must be unique, absent samples get ./.; 'records', like bcftools concat: identical sample lists, duplicate sites kept once), output (VCF path), include_records, overwrite.
    Merged QUAL = max, FILTER = PASS only if every contributing record passes (else the union of failing labels), INFO = first record carrying the key. Output is sorted.
    """

    options = dict(options or {})
    mode = _vcf.option_choice(options, 'mode', 'samples', ('samples', 'records'))
    merged = _vcf.merge_sets(_vcf.load_many(input), mode)
    return _done('merge_mt', input, options, _vcf_result(merged, options, mode=mode))

def merge_vcf_files(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Merge several VCF files, combining samples (default) or concatenating records.

    input: a list of at least two call sets (VCF text, VCF paths or {'records': [...]}) or {'inputs': [...]}.
    options: mode ('samples' default, like bcftools merge -m none: identical CHROM/POS/REF/ALT sites are combined, sample names must be unique, absent samples get ./.; 'records', like bcftools concat: identical sample lists, duplicate sites kept once), output (VCF path), include_records, overwrite.
    Merged QUAL = max, FILTER = PASS only if every contributing record passes (else the union of failing labels), INFO = first record carrying the key. Output is sorted.
    """

    options = dict(options or {})
    mode = _vcf.option_choice(options, 'mode', 'samples', ('samples', 'records'))
    merged = _vcf.merge_sets(_vcf.load_many(input), mode)
    return _done('merge_vcf_files', input, options, _vcf_result(merged, options, mode=mode))

def normalize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Normalise variant records: trim shared bases and, when a reference is given, left-align indels (bcftools norm style).

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: reference (plain sequence, FASTA text/path or {contig: sequence}), reference_start (1-based position of a plain sequence's first base, default 1), check_ref ('warn' default: records whose REF does not match are left unchanged; 'error'), split_multiallelic (False), sort (True), hgvs_accession (accession string or {contig: accession}) to add an HGVS g. string per ALT, output (VCF path), include_records, overwrite.
    Algorithm of Tan, Abecasis & Kang (2015): drop the shared last base, extend left by the reference base when an allele becomes empty, then drop shared leading bases. Multi-allelic records are handled jointly. Without options.reference only trimming (keeping one anchor base) is done and left-aligned=false is reported.
    """

    options = dict(options or {})
    return _done('normalize', input, options, _normalize_impl(input, options, require_reference=False))

def normalize_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Normalise a variant-by-sample (matrix-table style) call set: trim shared bases and, with a reference, left-align indels.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: reference (plain sequence, FASTA text/path or {contig: sequence}), reference_start (1-based position of a plain sequence's first base, default 1), check_ref ('warn' default: records whose REF does not match are left unchanged; 'error'), split_multiallelic (False), sort (True), hgvs_accession (accession string or {contig: accession}) to add an HGVS g. string per ALT, output (VCF path), include_records, overwrite.
    Algorithm of Tan, Abecasis & Kang (2015): drop the shared last base, extend left by the reference base when an allele becomes empty, then drop shared leading bases. Multi-allelic records are handled jointly. Without options.reference only trimming (keeping one anchor base) is done and left-aligned=false is reported.
    """

    options = dict(options or {})
    return _done('normalize_mt', input, options, _normalize_impl(input, options, require_reference=False))

def normalize_vcf(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Normalise a VCF: trim shared bases and, when a reference is given, left-align indels (bcftools norm style).

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: reference (plain sequence, FASTA text/path or {contig: sequence}), reference_start (1-based position of a plain sequence's first base, default 1), check_ref ('warn' default: records whose REF does not match are left unchanged; 'error'), split_multiallelic (False), sort (True), hgvs_accession (accession string or {contig: accession}) to add an HGVS g. string per ALT, output (VCF path), include_records, overwrite.
    Algorithm of Tan, Abecasis & Kang (2015): drop the shared last base, extend left by the reference base when an allele becomes empty, then drop shared leading bases. Multi-allelic records are handled jointly. Without options.reference only trimming (keeping one anchor base) is done and left-aligned=false is reported.
    """

    options = dict(options or {})
    return _done('normalize_vcf', input, options, _normalize_impl(input, options, require_reference=False))

def phase_genotypes(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Phase a child's heterozygous genotypes by transmission from the parents' genotypes (trio phasing).

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples). It must contain the child, father and mother samples.
    options: trio={'child','father','mother'} (or child/father/mother keys) sample names (required), output (optional VCF path), overwrite.
    The child GT becomes paternal|maternal when exactly one allele assignment is compatible with the parents; ambiguous, inconsistent and uninformative sites stay unphased. Statistical phasing (Beagle/SHAPEIT) is not performed.
    """

    options = dict(options or {})
    return _done('phase_genotypes', input, options, _analysis.phase_by_transmission(input, options))

def plot_density(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Plot variant density along chromosomes as a binned SVG bar chart and return the bin counts.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: bin_size (default 1000000), chroms (list), chrom_lengths ({contig: length}), svg_output (optional path to write the SVG), overwrite.
    The result has per-bin counts and density per Mb plus the SVG text.
    """

    options = dict(options or {})
    return _done('plot_density', input, options, _analysis.density_plot(input, options))

def plot_density_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Plot variant density of a variant-by-sample (matrix-table style) call set as a binned SVG bar chart.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: bin_size (default 1000000), chroms (list), chrom_lengths ({contig: length}), svg_output (optional path to write the SVG), overwrite.
    The result has per-bin counts and density per Mb plus the SVG text.
    """

    options = dict(options or {})
    return _done('plot_density_mt', input, options, _analysis.density_plot(input, options))

def prioritize_candidate_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Rank candidate variants with a transparent additive score from annotations already on the records (consequence, ClinVar, population AF, CADD/REVEL, gene panel).

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples). Annotations are read from record keys or INFO (consequence, clinvar/CLNSIG, gnomad_af, cadd, revel, gene).
    options: max_af (0.05; variants above it are excluded), exclude_benign (True), genes (panel list, +0.5), top_n.
    Heuristic ranking, not a validated classifier; every score component is returned. A variant with no population-frequency annotation keeps the full rarity score (0.5) but is flagged `frequency_unknown: true`, and one warning reports how many variants lacked frequency data (n_frequency_unknown).
    """

    options = dict(options or {})
    return _done('prioritize_candidate_variants', input, options, _prior.prioritize(input, options))

def split_multiallelic_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Split multi-allelic sites of a variant-by-sample (matrix-table style) call set into one record per ALT allele.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: output (VCF path), include_records, overwrite.
    GT is remapped (the selected ALT becomes 1, every other allele 0), FORMAT/AD and INFO AC/AF are subset to the selected allele and FORMAT/PL is dropped.
    """

    options = dict(options or {})
    vs = _vcf.load_variant_set(input)
    out = _vcf.split_multiallelic_set(vs)
    return _done('split_multiallelic_mt', input, options, _vcf_result(out, options, input_records=len(vs.records), split_records=len(out.records) - len(vs.records)))

def split_multiallelic_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Split multi-allelic variant records into one record per ALT allele.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: output (VCF path), include_records, overwrite.
    GT is remapped (the selected ALT becomes 1, every other allele 0), FORMAT/AD and INFO AC/AF are subset to the selected allele and FORMAT/PL is dropped.
    """

    options = dict(options or {})
    vs = _vcf.load_variant_set(input)
    out = _vcf.split_multiallelic_set(vs)
    return _done('split_multiallelic_variants', input, options, _vcf_result(out, options, input_records=len(vs.records), split_records=len(out.records) - len(vs.records)))

def summarize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise a call set: counts by variant type and chromosome, Ti/Tv, FILTER counts, multi-allelic sites and per-sample QC.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: none.
    Variants are classified per ALT allele by allele lengths after trimming (SNV, MNV, indel insertion/deletion/complex, SV for symbolic alleles or length change >= 50); Ti/Tv counts SNVs only; sample QC reports call rate, het/hom-alt counts and ratio, Ti/Tv and mean depth/GQ.
    """

    options = dict(options or {})
    return _done('summarize', input, options, _vcf.summarize_set(_vcf.load_variant_set(input)))

def summarize_mt(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Summarise a variant-by-sample (matrix-table style) call set: type counts, Ti/Tv, FILTER counts and per-sample QC.

    input: VCF text (tab-separated, optionally gzip bytes), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, optional qual/filter/info and genotypes/samples).
    options: none.
    Variants are classified per ALT allele by allele lengths after trimming (SNV, MNV, indel insertion/deletion/complex, SV for symbolic alleles or length change >= 50); Ti/Tv counts SNVs only; sample QC reports call rate, het/hom-alt counts and ratio, Ti/Tv and mean depth/GQ.
    """

    options = dict(options or {})
    return _done('summarize_mt', input, options, _vcf.summarize_set(_vcf.load_variant_set(input)))

__all__ = ['annotate', 'annotate_coding_effect', 'annotate_conservation_score', 'annotate_mt', 'annotate_pathogenicity_score', 'annotate_population_frequency', 'annotate_regulatory_effect', 'annotate_splice_effect', 'annotate_variant_consequence', 'annotate_vcf_ensembl_vep', 'annotate_vcf_vcfanno', 'calculate_allele_frequency', 'call_cutesv', 'call_deepvariant', 'call_delly', 'calculate_copy_number_profile', 'calculate_genotype_likelihood', 'calculate_haplotype_frequency', 'calculate_heterozygosity', 'calculate_hwe_pvalue', 'calculate_kinship_king', 'calculate_ld_matrix', 'calculate_minor_allele_frequency', 'calculate_missingness', 'calculate_runs_of_homozygosity', 'calculate_transition_transversion', 'calculate_variant_density', 'call', 'call_freebayes', 'call_gatk_haplotypecaller', 'call_sniffles2', 'call_mt', 'compare', 'compare_mt', 'compare_vcfs', 'construct_haplotypes', 'count', 'count_mt', 'detect_breakpoints', 'detect_low_frequency_variants', 'detect_mobile_element_insertions', 'detect_mosaic_variants', 'detect_repeat_expansions', 'filter', 'filter_mt', 'filter_variants_by_allele_frequency', 'filter_variants_by_depth', 'filter_variants_by_genotype_quality', 'filter_variants_by_missingness', 'filter_variants_by_quality', 'filter_vcf_vcftools', 'genotype_structural_variants', 'fill_missing_genotypes_with_mode', 'left_align', 'left_align_mt', 'left_align_variants', 'load_gemini_database', 'merge', 'merge_mt', 'merge_vcf_files', 'merge_structural_variants_survivor', 'normalize', 'normalize_mt', 'normalize_vcf', 'parse_vcf_records', 'phase_genotypes', 'plot_density', 'plot_density_mt', 'prioritize_candidate_variants', 'query_gemini_database', 'split_multiallelic', 'split_multiallelic_mt', 'split_multiallelic_variants', 'summarize', 'summarize_mt', 'read_vcf', 'fetch_vcf_region', 'calculate_observed_heterozygosity', 'calculate_fst', 'calculate_pca', 'calculate_tajimas_d', 'calculate_watterson_theta', 'calculate_nucleotide_diversity', 'calculate_site_frequency_spectrum', 'calculate_windowed_ld']

__candidate_adapter_operations__ = set(['annotate', 'annotate_coding_effect', 'annotate_conservation_score', 'annotate_mt', 'annotate_pathogenicity_score', 'annotate_population_frequency', 'annotate_regulatory_effect', 'annotate_splice_effect', 'annotate_variant_consequence', 'calculate_allele_frequency', 'calculate_copy_number_profile', 'calculate_genotype_likelihood', 'calculate_haplotype_frequency', 'calculate_heterozygosity', 'calculate_hwe_pvalue', 'calculate_kinship_king', 'calculate_ld_matrix', 'calculate_minor_allele_frequency', 'calculate_missingness', 'calculate_runs_of_homozygosity', 'calculate_transition_transversion', 'calculate_variant_density', 'call', 'call_mt', 'compare', 'compare_mt', 'compare_vcfs', 'construct_haplotypes', 'count', 'count_mt', 'detect_breakpoints', 'detect_low_frequency_variants', 'detect_mobile_element_insertions', 'detect_mosaic_variants', 'detect_repeat_expansions', 'filter', 'filter_mt', 'filter_variants_by_allele_frequency', 'filter_variants_by_depth', 'filter_variants_by_genotype_quality', 'filter_variants_by_missingness', 'filter_variants_by_quality', 'genotype_structural_variants', 'fill_missing_genotypes_with_mode', 'left_align', 'left_align_mt', 'left_align_variants', 'merge', 'merge_mt', 'merge_vcf_files', 'normalize', 'normalize_mt', 'normalize_vcf', 'parse_vcf_records', 'phase_genotypes', 'plot_density', 'plot_density_mt', 'prioritize_candidate_variants', 'split_multiallelic', 'split_multiallelic_mt', 'split_multiallelic_variants', 'summarize', 'summarize_mt'])

# These summaries are implemented locally and do not require a caller,
# annotation database, or optional backend.
__candidate_adapter_operations__ -= {
    'calculate_allele_frequency',
    'calculate_heterozygosity',
    'calculate_missingness',
    'calculate_transition_transversion',
    'calculate_variant_density',
}

# These record-level operations are deterministic and dependency-free. The
# broad wrappers above remain available for unimplemented operations, while
# these names resolve to the tested local implementations.
calculate_minor_allele_frequency = _calculate_minor_allele_frequency_local
compare = _compare_local
count = _count_local
filter_variants_by_allele_frequency = _filter_variants_by_allele_frequency_local
filter_variants_by_depth = _filter_variants_by_depth_local
filter_variants_by_genotype_quality = _filter_variants_by_genotype_quality_local
filter_variants_by_missingness = _filter_variants_by_missingness_local
filter_variants_by_quality = _filter_variants_by_quality_local
parse_vcf_records = _parse_vcf_records_local
split_multiallelic = _split_multiallelic_local
__candidate_adapter_operations__ -= {
    'calculate_minor_allele_frequency',
    'compare',
    'count',
    'filter_variants_by_allele_frequency',
    'filter_variants_by_depth',
    'filter_variants_by_genotype_quality',
    'filter_variants_by_missingness',
    'filter_variants_by_quality',
    'parse_vcf_records',
    'split_multiallelic',
}


# Optional-backend exports (see docs/api/features) 
from .nomenclature import (
    check_hgvs_syntax,
    describe_variant,
    format_hgvs_variant,
    load_vcf_variants,
    parse_hgvs_variant,
)

from .genotype_formats import (
    calculate_bgen_allele_frequency,
    read_bgen_probabilities,
    read_bgen_variants,
    read_plink_bed,
    read_plink_genotypes,
    read_plink_metadata,
    summarize_plink_bed,
    write_plink_bed,
    write_plink_genotypes,
)

__all__ = [*__all__, "check_hgvs_syntax", "describe_variant", "format_hgvs_variant", "load_vcf_variants", "parse_hgvs_variant", "calculate_bgen_allele_frequency", "read_bgen_probabilities", "read_bgen_variants", "read_plink_bed", "read_plink_genotypes", "read_plink_metadata", "summarize_plink_bed", "write_plink_bed", "write_plink_genotypes"]

# Optional-backend exports (hail)
from .hail_backends import (
    calculate_association_tests,
    calculate_hardy_weinberg_pvalues,
    calculate_hwe_normalized_pca,
    calculate_sample_qc,
    calculate_variant_qc,
    import_vcf_matrix_table,
    prune_linkage_disequilibrium,
    summarize_vcf_matrix_table,
)

__all__ = [*__all__, "calculate_association_tests", "calculate_hardy_weinberg_pvalues", "calculate_hwe_normalized_pca", "calculate_sample_qc", "calculate_variant_qc", "import_vcf_matrix_table", "prune_linkage_disequilibrium", "summarize_vcf_matrix_table"]

# Full-pass exports (see docs/api/features)
from .allel_functions import (
    calculate_ehh_decay,
    calculate_expected_heterozygosity,
    calculate_garud_h,
    calculate_hudson_fst_per_site,
    calculate_ihs,
    calculate_inbreeding_coefficient,
    calculate_mean_pairwise_difference,
    calculate_nsl,
    calculate_pairwise_distance,
    calculate_patterson_d,
    calculate_population_branching_statistic,
    calculate_rogers_huff_r,
    calculate_sequence_divergence,
    calculate_windowed_nucleotide_diversity,
    calculate_windowed_tajimas_d,
    calculate_xpehh,
    call_roh_hmm,
    find_mendel_errors,
    phase_trio_by_transmission,
    select_unlinked_variants,
)

from .bed_reader_functions import (
    calculate_genotype_relationship_matrix,
    standardize_genotypes,
    subset_genotype_matrix,
)

from .bgen_reader_functions import (
    create_bgen_metafile,
    read_bgen_allele_expectation,
    read_bgen_dosage,
)

from .cyvcf2_functions import (
    pearson_genotype_correlation,
    summarize_pairwise_relatedness,
    summarize_sample_heterozygosity,
)

from .hgvs_functions import (
    assemble_repeat_blocks,
    build_transcript_cigar,
    detect_repeat_blocks,
    normalize_indel_alleles,
    trim_common_prefix_strings,
    validate_hgvs_type_accession,
)

from .pandas_plink_functions import (
    read_gcta_relationship_matrix,
    read_plink_relationship_matrix,
)

from .pysam_functions import (
    annotate_vcf_from_tabix,
    build_vcf_consensus,
    compress_bgzip_file,
    concatenate_vcf_files,
    fetch_padded_sequence,
    index_tabix_file,
    intersect_vcf_sites,
    parse_genomic_region,
    parse_genomic_regions,
    query_vcf_fields,
    select_vcf_records,
    sort_vcf_records,
)

from .varcode_functions import (
    breakend_alt_sides,
    breakend_insert_bases,
    breakend_mate_position,
    check_vcf_sample_identity,
    classify_breakend_pair,
    compare_sv_calls,
    detect_germline_allele_overlap,
    enumerate_cryptic_exon_candidates,
    enumerate_phase_hypotheses,
    four_gamete_phase_test,
    load_maf_variants,
    normalize_nucleotide_allele,
    parse_symbolic_sv_alt,
    parse_vcf_genotype,
    partition_germline_by_phase,
    score_splice_site,
    trim_shared_allele_flanks,
    variant_overlaps_interval,
)

from .vcfpy_functions import (
    classify_vcf_alt,
    classify_vcf_substitution,
    format_vcf_field_value,
    parse_vcf_breakend,
    parse_vcf_field_value,
    remove_vcf_header_lines,
)

__all__ = [*__all__, "calculate_ehh_decay", "calculate_expected_heterozygosity", "calculate_garud_h", "calculate_hudson_fst_per_site", "calculate_ihs", "calculate_inbreeding_coefficient", "calculate_mean_pairwise_difference", "calculate_nsl", "calculate_pairwise_distance", "calculate_patterson_d", "calculate_population_branching_statistic", "calculate_rogers_huff_r", "calculate_sequence_divergence", "calculate_windowed_nucleotide_diversity", "calculate_windowed_tajimas_d", "calculate_xpehh", "call_roh_hmm", "find_mendel_errors", "phase_trio_by_transmission", "select_unlinked_variants", "calculate_genotype_relationship_matrix", "standardize_genotypes", "subset_genotype_matrix", "create_bgen_metafile", "read_bgen_allele_expectation", "read_bgen_dosage", "pearson_genotype_correlation", "summarize_pairwise_relatedness", "summarize_sample_heterozygosity", "assemble_repeat_blocks", "build_transcript_cigar", "detect_repeat_blocks", "normalize_indel_alleles", "trim_common_prefix_strings", "validate_hgvs_type_accession", "read_gcta_relationship_matrix", "read_plink_relationship_matrix", "annotate_vcf_from_tabix", "build_vcf_consensus", "compress_bgzip_file", "concatenate_vcf_files", "fetch_padded_sequence", "index_tabix_file", "intersect_vcf_sites", "parse_genomic_region", "parse_genomic_regions", "query_vcf_fields", "select_vcf_records", "sort_vcf_records", "breakend_alt_sides", "breakend_insert_bases", "breakend_mate_position", "check_vcf_sample_identity", "classify_breakend_pair", "compare_sv_calls", "detect_germline_allele_overlap", "enumerate_cryptic_exon_candidates", "enumerate_phase_hypotheses", "four_gamete_phase_test", "load_maf_variants", "normalize_nucleotide_allele", "parse_symbolic_sv_alt", "parse_vcf_genotype", "partition_germline_by_phase", "score_splice_site", "trim_shared_allele_flanks", "variant_overlaps_interval", "classify_vcf_alt", "classify_vcf_substitution", "format_vcf_field_value", "parse_vcf_breakend", "parse_vcf_field_value", "remove_vcf_header_lines"]

# Review-pass exports (see docs/api/features)
from .allel_review_functions import (
    build_equally_accessible_windows,
    calculate_between_population_pi,
    calculate_blockwise_patterson_f3,
    calculate_blockwise_patterson_fst,
    calculate_haplotype_diversity,
    calculate_moving_delta_tajimas_d,
    calculate_moving_haplotype_diversity,
    calculate_moving_patterson_f3,
    calculate_moving_patterson_fst,
    calculate_pairwise_dxy,
    calculate_patterson_f2,
    calculate_patterson_f3,
    calculate_patterson_fst,
    calculate_windowed_fixed_difference_density,
    calculate_windowed_patterson_fst,
    calculate_xpnsl,
    count_variants_in_windows,
    fold_joint_site_frequency_spectrum,
    locate_population_fixed_differences,
    locate_population_private_alleles,
    paint_parental_transmission,
    paint_shared_haplotype_prefixes,
    read_vcf_header_summary,
    scale_folded_joint_site_frequency_spectrum,
    scale_folded_site_frequency_spectrum,
    scale_joint_site_frequency_spectrum,
    scale_site_frequency_spectrum,
    standardize_score_by_allele_count,
    vcf_to_hdf5_archive,
    vcf_to_npz_archive,
    windowed_values_per_base,
)

from .hgvs_review_functions import (
    hgvs_interbase_span,
    repeat_notation_for_alleles,
    trim_common_suffix_strings,
)

from .pysam_review_functions import (
    add_bam_read_group,
    annotate_vcf_consequence,
    calculate_bam_read_depth,
    check_alignment_file_integrity,
    check_vcf_sample_concordance,
    clip_amplicon_primers,
    collate_bam_reads,
    concatenate_bam_files,
    count_bam_reads_per_interval,
    decode_sam_flags,
    detect_bcftools_roh,
    extract_fasta_from_bam,
    extract_fastq_from_bam,
    filter_vcf_by_expression,
    fix_bam_mate_fields,
    import_fastq_to_bam,
    index_alignment_file,
    index_fastq_file,
    index_vcf_csi,
    list_read_group_samples,
    mark_bam_duplicates,
    merge_bam_files,
    merge_indexed_vcfs,
    phase_bam_haplotypes,
    pileup_alignment_positions,
    pileup_bam_to_vcf,
    read_bam_header_and_records,
    read_vcf_header,
    recompute_nm_md_tags,
    reconstruct_reference_from_alignments,
    rename_vcf_samples,
    replace_bam_header,
    reset_alignments_to_unmapped,
    select_alignment_records,
    sort_alignment_file,
    split_bam_by_read_group,
    summarize_alignment_stats,
    summarize_amplicon_stats,
    summarize_bam_coverage,
    summarize_bam_index_stats,
    write_fasta_sequence_dictionary,
)

from .varcode_review_functions import (
    breakend_junction_span,
    has_only_standard_nucleotides,
    infer_reference_assembly,
    is_mitochondrial_contig,
    is_purine_base,
    is_ucsc_reference_name,
    join_breakend_layouts,
    pair_breakend_mates,
    project_intervals_onto_layout,
    rewrite_vcf_records,
)

__all__ = [*__all__, "build_equally_accessible_windows", "calculate_between_population_pi", "calculate_blockwise_patterson_f3", "calculate_blockwise_patterson_fst", "calculate_haplotype_diversity", "calculate_moving_delta_tajimas_d", "calculate_moving_haplotype_diversity", "calculate_moving_patterson_f3", "calculate_moving_patterson_fst", "calculate_pairwise_dxy", "calculate_patterson_f2", "calculate_patterson_f3", "calculate_patterson_fst", "calculate_windowed_fixed_difference_density", "calculate_windowed_patterson_fst", "calculate_xpnsl", "count_variants_in_windows", "fold_joint_site_frequency_spectrum", "locate_population_fixed_differences", "locate_population_private_alleles", "paint_parental_transmission", "paint_shared_haplotype_prefixes", "read_vcf_header_summary", "scale_folded_joint_site_frequency_spectrum", "scale_folded_site_frequency_spectrum", "scale_joint_site_frequency_spectrum", "scale_site_frequency_spectrum", "standardize_score_by_allele_count", "vcf_to_hdf5_archive", "vcf_to_npz_archive", "windowed_values_per_base", "hgvs_interbase_span", "repeat_notation_for_alleles", "trim_common_suffix_strings", "add_bam_read_group", "annotate_vcf_consequence", "calculate_bam_read_depth", "check_alignment_file_integrity", "check_vcf_sample_concordance", "clip_amplicon_primers", "collate_bam_reads", "concatenate_bam_files", "count_bam_reads_per_interval", "decode_sam_flags", "detect_bcftools_roh", "extract_fasta_from_bam", "extract_fastq_from_bam", "filter_vcf_by_expression", "fix_bam_mate_fields", "import_fastq_to_bam", "index_alignment_file", "index_fastq_file", "index_vcf_csi", "list_read_group_samples", "mark_bam_duplicates", "merge_bam_files", "merge_indexed_vcfs", "phase_bam_haplotypes", "pileup_alignment_positions", "pileup_bam_to_vcf", "read_bam_header_and_records", "read_vcf_header", "recompute_nm_md_tags", "reconstruct_reference_from_alignments", "rename_vcf_samples", "replace_bam_header", "reset_alignments_to_unmapped", "select_alignment_records", "sort_alignment_file", "split_bam_by_read_group", "summarize_alignment_stats", "summarize_amplicon_stats", "summarize_bam_coverage", "summarize_bam_index_stats", "write_fasta_sequence_dictionary", "breakend_junction_span", "has_only_standard_nucleotides", "infer_reference_assembly", "is_mitochondrial_contig", "is_purine_base", "is_ucsc_reference_name", "join_breakend_layouts", "pair_breakend_mates", "project_intervals_onto_layout", "rewrite_vcf_records"]

# Family-pass exports (see docs/api/features)
from .allel_family_functions import (
    calculate_blockwise_patterson_d,
    calculate_joint_site_frequency_spectrum_from_counts,
    calculate_moving_garud_h,
    calculate_moving_hudson_fst,
    calculate_moving_patterson_d,
    calculate_moving_tajimas_d,
    calculate_moving_weir_cockerham_fst,
    calculate_randomized_pca,
    calculate_rogers_huff_r_between,
    calculate_site_frequency_spectrum_from_counts,
    calculate_weir_cockerham_fst_per_variant,
    calculate_windowed_divergence,
    calculate_windowed_hudson_fst,
    calculate_windowed_watterson_theta,
    calculate_windowed_weir_cockerham_fst,
    phase_trio_step_by_transmission,
    read_vcf_field_table,
    vcf_to_csv_table,
)

from .pysam_family_functions import (
    export_vcf_to_phasing_format,
)

from .varcode_family_functions import (
    load_maf_table,
    read_vcf_table,
)

__all__ = [*__all__, "calculate_blockwise_patterson_d", "calculate_joint_site_frequency_spectrum_from_counts", "calculate_moving_garud_h", "calculate_moving_hudson_fst", "calculate_moving_patterson_d", "calculate_moving_tajimas_d", "calculate_moving_weir_cockerham_fst", "calculate_randomized_pca", "calculate_rogers_huff_r_between", "calculate_site_frequency_spectrum_from_counts", "calculate_weir_cockerham_fst_per_variant", "calculate_windowed_divergence", "calculate_windowed_hudson_fst", "calculate_windowed_watterson_theta", "calculate_windowed_weir_cockerham_fst", "phase_trio_step_by_transmission", "read_vcf_field_table", "vcf_to_csv_table", "export_vcf_to_phasing_format", "load_maf_table", "read_vcf_table"]

_register_public_verbose(globals())
