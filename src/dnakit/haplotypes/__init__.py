"""Haplotype operations: native trio/pedigree logic, haplotype statistics and LD tools.

Operations that need a reference panel or a statistical-genetics backend return
``status: "not_implemented"`` (see ``dnakit._candidate_runtime``).
"""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import (
    candidate_not_implemented,
    candidate_ok,
    run_candidate_operation,
)

from . import _impl_haplotypes as _impl

_MODULE = "haplotypes"
_TRIO_FORM = "input: {'child': [...], 'father': [...], 'mother': [...]} (+ optional 'sites' labels), one call per site"
_CALLS = "A call is 'a/b' or 'a|b' (integer alleles), a pair (a, b), an alt dosage 0/1/2, or None / './.' for missing."


def _ok(
    operation: str,
    input: object,
    result: object,
    options: Mapping[str, object] | None,
    warnings: tuple[str, ...] = (),
) -> object:
    return candidate_ok(_MODULE, operation, input, result, options=options, warnings=warnings)


def analyze_trio_inheritance(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Check Mendelian consistency of a child's genotypes against both parents and report the transmitted alleles.

    input: {'child': [...], 'father': [...], 'mother': [...]} with one genotype call per site (optional
    'sites' labels). A call is 'a/b' or 'a|b' (integer alleles), a pair (a, b), an alt dosage 0/1/2,
    or None / './.' for missing. A site is consistent if the child's two alleles can be formed by one
    allele from each parent; sites with a missing call are skipped. Returns per-site status
    (consistent / inconsistent / missing), the unique (paternal, maternal) allele assignment when
    determinable, counts and the Mendelian error rate = inconsistent / informative sites. options: none.
    """

    return _ok("analyze_trio_inheritance", input, _impl.trio_inheritance(input), options)


def calculate_pedigree_likelihood(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Compute the likelihood of observed biallelic genotypes on a small pedigree (exact enumeration, no genotyping error).

    input: {'pedigree': [{'id', 'father', 'mother'}, ...], 'genotypes': {id: [calls per site]}}; give both
    parents or neither (founders); members missing from 'genotypes' or with missing calls are summed out
    (at most 8 ungenotyped members per site). Founders follow Hardy-Weinberg with the alt-allele frequency
    p; transmission is Mendelian; sites are independent. Calls are 0/1 alleles ('0/1') or dosages 0/1/2.
    options: allele_frequency (number in (0,1) or one per site, default 0.5). Returns per-site and total
    likelihood and log10 likelihood (None if the data are impossible, likelihood 0).
    """

    settings = _impl.options_dict(options)
    return _ok(
        "calculate_pedigree_likelihood", input, _impl.pedigree_likelihood(input, settings), options
    )


def cluster(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Cluster haplotypes by single linkage on Hamming distance (connected components within max_distance mismatches).

    input: list of haplotype strings (equal length) or {name: haplotype | [haplotypes]}.
    options: max_distance (int >= 0, default 1); ignore_missing (default True: '?', 'N', '.', '-' are skipped
    when comparing a pair). At most 5,000 haplotypes. Returns clusters in order of their first member.
    """

    settings = _impl.options_dict(options)
    distance = _impl.get_int(settings, "max_distance", 1, low=0, high=10**6)
    return _ok(
        "cluster",
        input,
        _impl.cluster_haplotypes(input, distance or 0, bool(settings.get("ignore_missing", True))),
        options,
    )


def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare haplotypes by pairwise Hamming distance, or two phased haplotype pairs by switch error rate.

    options.mode 'hamming' (default): input is a list of at least two haplotype strings (equal length) or
    {name: haplotype | [haplotypes]}; returns the distance matrix (mismatching sites; '?', 'N', '.', '-'
    skipped unless ignore_missing=False) and, for exactly two haplotypes, the differing positions.
    options.mode 'switch_error': input {'truth': [h1, h2], 'inferred': [h1, h2]}; among sites heterozygous
    and identical as allele sets in both, a switch is a change of relative phase between consecutive sites and
    switch_error_rate = switches / (sites - 1).
    """

    settings = _impl.options_dict(options)
    return _ok("compare", input, _impl.compare_haplotypes(input, settings), options)


def construct(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Construct haplotype blocks from phased 0/1 haplotypes by a complete pairwise-LD rule.

    input: list of haplotypes (rows = chromosomes, equal-length strings or lists of 0/1 alleles) or
    {name: haplotype(s)}. Rule: left to right, a block is the longest run of at least min_snps polymorphic
    sites in which EVERY pair has LD (|D'| or r^2) >= threshold, LD computed from the given chromosomes
    (not Gabriel's confidence-interval method; small samples inflate D'). Monomorphic sites (and sites with
    MAF < min_maf) are excluded. options: statistic ('dprime' default or 'r2'); threshold (default 0.8);
    min_snps (default 2); min_maf (default 0.0); positions (optional coordinate per site, echoed for blocks).
    """

    settings = _impl.options_dict(options)
    return _ok("construct", input, _impl.ld_blocks(input, settings), options)


def count(input: object, *, options: Mapping[str, object] | None = None) -> object:
    return run_candidate_operation("haplotypes", "count", input, options=options)


def detect_compound_heterozygotes(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Find genes carrying heterozygous variants on both parental haplotypes (compound heterozygotes, in trans).

    input: {'genes': [gene per site], 'child': [...]} plus either 'father' and 'mother' (trio mode: the
    origin of a heterozygous child allele is taken from the unique Mendelian assignment of child alleles to
    parents) or no parents with PHASED child calls like '0|1' (phase mode: haplotype 1 vs 2). Calls: 'a/b',
    'a|b', pair, dosage, None. A gene is reported when it has at least min_variants_per_haplotype
    non-reference variants on each haplotype. Sites with ambiguous origin are listed as undetermined.
    options: min_variants_per_haplotype (int, default 1).
    """

    settings = _impl.options_dict(options)
    return _ok(
        "detect_compound_heterozygotes", input, _impl.compound_hets(input, settings), options
    )


def detect_denovo_variants(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Detect candidate de novo variants in a trio: child alleles absent from both parents' genotypes.

    input: {'child': [...], 'father': [...], 'mother': [...]} (+ optional 'sites' labels), one call per site.
    A site is a de novo candidate when the child carries an allele that neither parent's genotype contains
    (for example parents 0/0 and child 0/1); other Mendelian violations are listed separately. Sites with
    any missing call are skipped. Depth/quality filtering is not applied. options: none.
    """

    return _ok("detect_denovo_variants", input, _impl.denovo_sites(input), options)


def frequency(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Count haplotypes and compute their frequencies and Nei's haplotype diversity.

    input: list of haplotype strings, {sample: [h1, h2]} / {name: haplotype}, or {haplotype: count}.
    frequency = count / total chromosomes; haplotype diversity = n/(n-1) * (1 - sum p_i^2) (None if n < 2).
    Rows are sorted by decreasing count, then haplotype string. options: none.
    """

    return _ok("frequency", input, _impl.frequency_table(input), options)


def impute(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Not implemented: genotype imputation needs a phased reference panel and an HMM imputation engine."""

    return candidate_not_implemented(
        _MODULE,
        "impute",
        input,
        options=options,
        requires="phased reference haplotype panel and an imputation backend (Beagle, IMPUTE5, Minimac4)",
        message="Imputation infers missing genotypes from a reference panel with a statistical model; neither is part of the input.",
    )


def infer_parentage(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Screen candidate parents for a child by Mendelian exclusion (sites where parent and child share no allele).

    input: {'child': [calls], 'candidates': {name: [calls]}}, one call per site (see analyze_trio_inheritance for
    call formats). A candidate is excluded if the number of compared sites with no shared allele (for
    example opposite homozygotes) exceeds max_mismatches; the best candidates are those not excluded with the
    fewest mismatches. A candidate with no jointly called site (compared_sites 0) cannot be assessed: it is neither
    excluded nor in not_excluded/best_candidates but is listed in not_assessed, with a warning.
    This is exclusion only, not a paternity likelihood. options: max_mismatches (int, default 0).
    """

    settings = _impl.options_dict(options)
    result = _impl.infer_parentage(input, settings)
    warnings = tuple(result.pop("warnings", ()))
    return _ok("infer_parentage", input, result, options, warnings)


def phase(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Phase a child's genotypes using Mendelian trio logic (paternal|maternal order).

    input: {'child': [...], 'father': [...], 'mother': [...]}, one call per site. A heterozygous child site is
    phased when exactly one assignment of its two alleles to (father, mother) is compatible with the parents'
    genotypes; sites where both parents are heterozygous stay unphased ('a/b', status 'ambiguous'), Mendelian
    violations are flagged. Statistical (population) phasing and read-based phasing are not performed. options: none.
    """

    return _ok("phase", input, _impl.trio_phase(input), options)


def reconstruct(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate haplotype frequencies and each individual's haplotype pair from unphased biallelic genotypes by EM.

    input: {sample: [calls per site]} with calls as dosage 0/1/2 or '0/1' strings, no missing calls, at most 12
    heterozygous sites per individual. Excoffier-Slatkin (1995) EM under Hardy-Weinberg, uniform start;
    each individual gets its maximum-posterior haplotype pair. The reported log_likelihood is the log-likelihood
    of the data at the FINAL (converged) haplotype frequencies. options: max_iterations (default 1000);
    tolerance (max frequency change to stop, default 1e-8).
    """

    settings = _impl.options_dict(options)
    return _ok("reconstruct", input, _impl.reconstruct_em(input, settings), options)


def select(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Select haplotypes by minimum frequency, minimum count and/or top-N from a haplotype sample.

    input: as for frequency (list of haplotypes, {sample: [h1, h2]} or {haplotype: count}).
    options: min_frequency (default 0.05); min_count (default 1); top_n (optional). Returns the retained rows
    sorted by decreasing frequency and the total frequency they cover.
    """

    settings = _impl.options_dict(options)
    return _ok("select", input, _impl.select_haplotypes(input, settings), options)


def simulate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Simulate diploid individuals by drawing two haplotypes from given haplotype frequencies (random mating).

    input: {haplotype: frequency or count}. options: n_individuals (default 10, max 100000); seed (default 0;
    results are reproducible for a fixed seed). No mutation, recombination or demography is simulated;
    this is not a coalescent simulator.
    """

    settings = _impl.options_dict(options)
    return _ok("simulate", input, _impl.simulate_haplotypes(input, settings), options)


def tag(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Select tag SNPs from phased 0/1 haplotypes by greedy r^2 binning (Carlson et al. 2004).

    input: list of haplotypes (equal-length 0/1 strings or lists) or {name: haplotype(s)}. Repeatedly
    chooses the site whose r^2 with the most still-untagged sites is >= r2_threshold (ties: lowest index)
    and assigns those sites to its bin. Monomorphic sites are excluded. options: r2_threshold (default 0.8);
    min_maf (default 0.0). Returns tag site indices and their bins.
    """

    settings = _impl.options_dict(options)
    return _ok("tag", input, _impl.tag_snps(input, settings), options)


def visualize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw haplotypes as an SVG allele grid (rows = haplotypes, columns = sites; allele 1 dark, 0 light, other grey).

    input: list of haplotype strings (equal length) or {name: haplotype(s)}. options: output (optional file
    path; the SVG text is written there). Result holds the SVG string.
    """

    settings = _impl.options_dict(options)
    svg, rows, sites = _impl.haplotype_svg(input)
    return _ok(
        "visualize",
        input,
        {
            "svg": svg,
            "haplotype_count": rows,
            "site_count": sites,
            "output_path": _impl.write_output(settings, svg),
        },
        options,
    )


__all__ = [
    "analyze_trio_inheritance",
    "calculate_pedigree_likelihood",
    "cluster",
    "compare",
    "construct",
    "count",
    "detect_compound_heterozygotes",
    "detect_denovo_variants",
    "frequency",
    "impute",
    "infer_parentage",
    "phase",
    "reconstruct",
    "select",
    "simulate",
    "tag",
    "visualize",
]

__candidate_adapter_operations__ = set(
    [
        "analyze_trio_inheritance",
        "calculate_pedigree_likelihood",
        "cluster",
        "compare",
        "construct",
        "count",
        "detect_compound_heterozygotes",
        "detect_denovo_variants",
        "frequency",
        "impute",
        "infer_parentage",
        "phase",
        "reconstruct",
        "select",
        "simulate",
        "tag",
        "visualize",
    ]
)


from dnakit._verbosity import register_public as _register_public_verbose

_register_public_verbose(globals())
