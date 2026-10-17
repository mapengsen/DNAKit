"""Local candidate adapters for the clinical capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import candidate_not_implemented, candidate_ok
from dnakit.variants import _impl_prioritize as _prior

from . import _impl_clinical as _clin


def _done(operation: str, input: object, options: Mapping[str, object], result: dict[str, object]) -> dict[str, object]:
    """Wrap a computed result in the candidate envelope (``warnings`` move to the envelope)."""

    warnings = list(result.pop("warnings", []) or [])
    return candidate_ok("clinical", operation, input, result, options=options, warnings=warnings)


def analyze_trio_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Classify the variants of a child-father-mother trio by inheritance pattern: de novo, Mendelian error, inherited from one or both parents, homozygous with carrier parents, and compound-heterozygous candidates.

    input: VCF text (tab-separated), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, genotypes/samples with FORMAT fields). It needs the three trio samples.
    options: child, father, mother (sample names; if absent and the input has exactly three samples the order child, father, mother is assumed), min_depth (0), min_genotype_quality (0), require_quality_pass (True). Gene symbols are read from a record 'gene' key or INFO GENE/SYMBOL/GENEINFO for compound heterozygosity.
    Multi-allelic sites are split first. Autosomal diploid logic only; X-linked and mitochondrial inheritance are not modelled.
    """

    options = dict(options or {})
    return _done('analyze_trio_variants', input, options, _clin.analyze_trio(input, options))

def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate clinical variants with ClinVar, gnomAD and gene-disease information (not implemented).

    Needs the annotation databases; use dnakit.databases / dnakit.annotation look-ups for individual variants.

    Currently returns status 'not_implemented': clinical annotation needs reference databases that are not bundled.
    """

    return candidate_not_implemented('clinical', 'annotate', input, options=options,
        requires='ClinVar / gnomAD / gene-disease databases (local files or network access)',
        message='clinical annotation needs reference databases that are not bundled.')

def build_cancer_evolution_tree(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a rooted perfect-phylogeny evolution tree of somatic mutations from a binary mutation-by-sample (or cell) matrix.

    input: {'mutations': {'TP53': [1,1,0], ...}, 'samples': [...]} (or mutation -> {sample: 0/1}).
    options: none.
    Mutations with identical sample sets share a node; a node's parent is the node with the smallest strictly larger sample set; if two mutations' sample sets partially overlap (four-gamete violation) the result reports perfect_phylogeny=false with the violations and no tree. The infinite-sites assumption and error-free calls are required.
    """

    options = dict(options or {})
    return _done('build_cancer_evolution_tree', input, options, _clin.evolution_tree(input, options))

def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate an unspecified clinical quantity (not implemented).

    The name does not identify a computation; use a specific function (classify, qc, detect_de_novo_variants, ...).

    Currently returns status 'not_implemented': 'calculate' is too unspecific to have a defined meaning; call a specific clinical function.
    """

    return candidate_not_implemented('clinical', 'calculate', input, options=options,
        requires='a specific clinical calculation (this generic name has no defined meaning)',
        message="'calculate' is too unspecific to have a defined meaning; call a specific clinical function.")

def call(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Call clinical variants from sequencing reads (not implemented).

    Needs a variant caller; use dnakit.variants.call_freebayes / call_gatk_haplotypecaller / call_deepvariant.

    Currently returns status 'not_implemented': variant calling needs an external caller and aligned reads.
    """

    return candidate_not_implemented('clinical', 'call', input, options=options,
        requires='variant caller binary and aligned reads',
        message='variant calling needs an external caller and aligned reads.')

def classify(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Combine ACMG/AMP 2015 evidence criteria into a five-tier variant classification (Pathogenic ... Benign).

    input: the criteria met, e.g. ['PVS1', 'PM2_Supporting'] or {'criteria': [...], 'id': ...} or {'variants': [{'id', 'criteria'}, ...]}; entries may be codes or {'code', 'strength'}.
    options: none.
    Applies Table 5 of Richards et al. 2015: Pathogenic (PVS1+>=1 PS, PVS1+>=2 PM, PVS1+PM+PP, PVS1+>=2 PP, >=2 PS, PS+>=3 PM, PS+2 PM+>=2 PP, PS+PM+>=4 PP); Likely pathogenic (PVS1+PM, PS+1-2 PM, PS+>=2 PP, >=3 PM, 2 PM+>=2 PP, PM+>=4 PP); Benign (BA1 or >=2 BS); Likely benign (BS+BP or >=2 BP); otherwise, or when pathogenic and benign rules are both met, Uncertain significance. Codes: PVS1, PS1-4, PM1-6, PP1-5, BA1, BS1-4, BP1-7; a suffix (PM2_Supporting) changes the strength. The criteria themselves are not evaluated from data.
    """

    options = dict(options or {})
    return _done('classify', input, options, _clin.classify_acmg(input, options))

def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare clinical data sets (not implemented).

    The name does not identify a clinical comparison; for call-set comparison use dnakit.variants.compare_vcfs.

    Currently returns status 'not_implemented': 'compare' is too unspecific here; use dnakit.variants.compare_vcfs for variant call sets.
    """

    return candidate_not_implemented('clinical', 'compare', input, options=options,
        requires='a specific comparison definition (this generic name has no defined meaning)',
        message="'compare' is too unspecific here; use dnakit.variants.compare_vcfs for variant call sets.")

def detect_circulating_tumor_dna(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Test for circulating tumour DNA in plasma at tracked tumour variants (tumour-informed detection by a pooled binomial test).

    input: tracked sites [{'alt_count','depth'[,'id','error_rate']}] or {'sites': [...]}, or a VCF/records with FORMAT/AD for the plasma sample (options.sample).
    options: background_error_rate (0.001, must be measured for the assay), alpha (0.01), site_alpha (0.01), min_alt_reads (2), min_positive_sites (1), vaf_to_fraction_factor (2.0).
    Detected when the pooled ALT count is improbable under the background error rate (one-sided binomial p <= alpha) and enough sites are individually positive; the tumour fraction is estimated as factor*max(0, pooled VAF - error). Not a validated clinical assay.
    """

    options = dict(options or {})
    return _done('detect_circulating_tumor_dna', input, options, _clin.tumor_informed_detection(input, options, alpha=0.01, min_sites=1, label='ctDNA detection'))

def detect_de_novo_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Detect de novo variant candidates in a child-father-mother trio from genotypes, depths and allelic depths.

    input: VCF text (tab-separated), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, genotypes/samples with FORMAT fields). It needs the three trio samples.
    options: child, father, mother (sample names), min_child_depth (10), min_parent_depth (10), min_genotype_quality (20), min_child_vaf (0.25, applied when AD exists), max_parent_alt_reads (0, applied when AD exists).
    A candidate is a site where the child carries an ALT allele and both parents are called homozygous without it; candidates failing the filters are counted in rejected_by, once per variant, under the FIRST failing check in a fixed order: child/father/mother depth, child/father/mother GQ, child ALT fraction (min_child_vaf), father/mother ALT reads (max_parent_alt_reads). The genotype checks (child carries ALT, both parents called homozygous without it) come first and define the candidates.
    """

    options = dict(options or {})
    return _done('detect_de_novo_variants', input, options, _clin.de_novo(input, options))

def detect_minimal_residual_disease(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Test for minimal residual disease (MRD) from plasma read counts at tracked tumour variants with a stricter tumour-informed pooled binomial test.

    input: tracked sites [{'alt_count','depth'[,'id','error_rate']}] or {'sites': [...]}, or a VCF/records with FORMAT/AD for the plasma sample (options.sample).
    options: background_error_rate (0.001), alpha (0.001), site_alpha (0.01), min_alt_reads (2), min_positive_sites (2), vaf_to_fraction_factor (2.0).
    Same statistics as detect_circulating_tumor_dna with stricter defaults (alpha 0.001 and at least two positive sites). Not a validated clinical assay.
    """

    options = dict(options or {})
    return _done('detect_minimal_residual_disease', input, options, _clin.tumor_informed_detection(input, options, alpha=0.001, min_sites=2, label='MRD detection'))

def generate_clinical_variant_report(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Generate a clinical variant report (Markdown or plain text) from classified variants.

    input: {'variants': [{'gene','hgvs','zygosity','classification' or 'criteria', 'clinvar', ...}], 'sample': ..., 'indication': ...}, a list of such rows, or a VCF/records (gene, ClinVar and consequence fields are read from the records).
    options: title, format ('markdown' default or 'text'), include_vus (True), hgvs_accession (string or {contig: accession}) and reference to build HGVS g. strings for VCF input, output (path), overwrite.
    Classifications may be given as Pathogenic/Likely pathogenic/Uncertain significance/Likely benign/Benign (or P/LP/VUS/LB/B); when only 'criteria' is given they are combined with the ACMG/AMP 2015 rules. The report lists reportable findings, VUS and counts; it is a summary of the supplied evidence, not a validated diagnostic report.
    """

    options = dict(options or {})
    return _done('generate_clinical_variant_report', input, options, _clin.clinical_report(input, options))

def infer_tumor_clones(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Infer tumour clones by clustering mutation ALT fractions with a binomial mixture model (EM, number of clones chosen by BIC).

    input: mutations [{'alt_count','depth'[,'id','copy_number','multiplicity']}] or {'sites': [...]}, or a VCF/records with FORMAT/AD (options.sample).
    options: purity (1.0), max_clones (4, max 8), copy_number (default 2), sample.
    CCF = VAF*(purity*CN + 2*(1-purity))/(purity*multiplicity) per mutation; clusters are ranked by ALT fraction (the highest is called clonal). Copy-number changes and clone phylogeny are not modelled.
    """

    options = dict(options or {})
    return _done('infer_tumor_clones', input, options, _clin.infer_clones(input, options))

def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Predict clinical outcomes or variant effects with a trained model (not implemented).

    Needs a trained predictor and its weights.

    Currently returns status 'not_implemented': prediction needs trained model weights that are not available.
    """

    return candidate_not_implemented('clinical', 'predict', input, options=options,
        requires='trained clinical/pathogenicity prediction model',
        message='prediction needs trained model weights that are not available.')

def prioritize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Rank clinical candidate variants with a transparent additive score (consequence, ClinVar, rarity, CADD/REVEL, gene panel).

    input: VCF text (tab-separated), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, genotypes/samples with FORMAT fields). Annotations are read from record keys or INFO (consequence, clinvar/CLNSIG, gnomad_af, cadd, revel, gene).
    options: max_af (0.05; variants above it are excluded), exclude_benign (True), genes (panel list, +0.5), top_n.
    Heuristic ranking, not a validated classifier; every score component is returned. A variant with no population-frequency annotation keeps the full rarity score (0.5) but is flagged `frequency_unknown: true`, and one warning reports how many variants lacked frequency data (n_frequency_unknown).
    """

    options = dict(options or {})
    return _done('prioritize', input, options, _prior.prioritize(input, options))

def qc(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Quality-control a clinical call set: per-sample call rate, het/hom ratio, Ti/Tv and depth, per-site call rate and exact HWE test, with threshold flags.

    input: VCF text (tab-separated), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, genotypes/samples with FORMAT fields).
    options: min_sample_call_rate (0.95), min_site_call_rate (0.9), min_hwe_pvalue (1e-6), min_mean_depth, max_het_hom_alt_ratio, ti_tv_range [low, high] (the last three are off unless given).
    This QC covers called genotypes only (not raw sequencing-run metrics).
    """

    options = dict(options or {})
    return _done('qc', input, options, _clin.clinical_qc(input, options))

def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Read a (ClinVar-style) VCF or variant records into clinical fields: gene, ClinVar significance and review status, disease, consequence and genotypes.

    input: VCF text (tab-separated), a path to a VCF/VCF.gz, or a list of variant mappings (chrom, pos, ref, alt, genotypes/samples with FORMAT fields).
    options: none.
    Reads INFO CLNSIG, CLNREVSTAT, CLNDN, GENEINFO/GENE and MC; significance strings are also normalised into lower-case terms. No external look-up is made.
    """

    options = dict(options or {})
    return _done('read', input, options, _clin.read_clinical(input, options))

def report(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Generate a clinical variant report (Markdown or plain text) from classified variants.

    input: {'variants': [{'gene','hgvs','zygosity','classification' or 'criteria', 'clinvar', ...}], 'sample': ..., 'indication': ...}, a list of such rows, or a VCF/records (gene, ClinVar and consequence fields are read from the records).
    options: title, format ('markdown' default or 'text'), include_vus (True), hgvs_accession (string or {contig: accession}) and reference to build HGVS g. strings for VCF input, output (path), overwrite.
    Classifications may be given as Pathogenic/Likely pathogenic/Uncertain significance/Likely benign/Benign (or P/LP/VUS/LB/B); when only 'criteria' is given they are combined with the ACMG/AMP 2015 rules. The report lists reportable findings, VUS and counts; it is a summary of the supplied evidence, not a validated diagnostic report.
    """

    options = dict(options or {})
    return _done('report', input, options, _clin.clinical_report(input, options))

def type_hla(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Type HLA alleles from sequencing reads (not implemented).

    Needs an HLA typing tool (OptiType, HLA-LA, arcasHLA, T1K) plus the IPD-IMGT/HLA reference.

    Currently returns status 'not_implemented': HLA typing needs an external program and the IMGT/HLA reference.
    """

    return candidate_not_implemented('clinical', 'type_hla', input, options=options,
        requires='HLA typing program (OptiType/HLA-LA/arcasHLA) with IPD-IMGT/HLA reference and reads',
        message='HLA typing needs an external program and the IMGT/HLA reference.')

__all__ = ['analyze_trio_variants', 'annotate', 'build_cancer_evolution_tree', 'calculate', 'call', 'classify', 'compare', 'detect_circulating_tumor_dna', 'detect_de_novo_variants', 'detect_minimal_residual_disease', 'generate_clinical_variant_report', 'infer_tumor_clones', 'predict', 'prioritize', 'qc', 'read', 'report', 'type_hla']

__candidate_adapter_operations__ = set(['analyze_trio_variants', 'annotate', 'build_cancer_evolution_tree', 'calculate', 'call', 'classify', 'compare', 'detect_circulating_tumor_dna', 'detect_de_novo_variants', 'detect_minimal_residual_disease', 'generate_clinical_variant_report', 'infer_tumor_clones', 'predict', 'prioritize', 'qc', 'read', 'report', 'type_hla'])


from dnakit._verbosity import register_public as _register_public_verbose
_register_public_verbose(globals())
