"""Forensic and reproductive-genetics statistics from allele frequencies and single-source STR/SNP
profiles.

Profiles are ``{locus: [allele1, allele2]}``; allele frequencies are ``{locus: {allele:
frequency}}`` (or records with columns locus, allele, frequency).  Formulas are stated in each
result.  Every function returns the DNAKit candidate envelope.
"""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import candidate_not_implemented
from dnakit.gwas._impl_common import finish as _finish
from dnakit.gwas._impl_common import get_options as _get_options

from . import _impl_forensic as _f

_M = "forensic_reproductive"


def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate a forensic statistic from allele frequencies: random match probability,
    likelihood ratio, paternity index or kinship index.

    options: quantity = 'random_match_probability' (default; input {'profile',
    'frequencies'}: product rule p^2 / 2pq per locus, optional theta with the NRC II 4.10
    correction), 'likelihood_ratio' (input {'evidence', 'suspect', 'frequencies'}; LR =
    1/P(genotype) per matching locus, 0 on a mismatch), 'paternity_index' (input {'child',
    'alleged_father', 'mother' (optional), 'frequencies'}; PI = P(child | mother,
    father)/P(child | mother, random man), combined index and probability of paternity for
    prior=0.5; no mutation model) or 'kinship_index' (input {'profile_1', 'profile_2',
    'frequencies'}; relationship in parent_child, full_sibling, half_sibling,
    grandparent_grandchild, uncle_nephew, first_cousin, identical_twin or custom
    ibd=[k0,k1,k2]). Alleles absent from the frequency table raise an error unless
    min_frequency is set.
    """

    return _finish(_M, "calculate", input, options, _f.calculate)


def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare two single-source genotype profiles locus by locus (identical / one allele shared
    / no allele shared).

    input: {'profile_1': {...}, 'profile_2': {...}} and optionally 'frequencies' (adds the
    random match probability of profile_1). Reports per-locus status, counts, loci present in
    only one profile and whether all compared loci are identical.
    """

    return _finish(_M, "compare", input, options, _f.compare)


def deconvolute(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Deconvolute a DNA mixture into contributor genotypes (not implemented: needs a
    probabilistic genotyping program).
    """

    return candidate_not_implemented(
        _M,
        "deconvolute",
        input,
        options=_get_options(options),
        requires=(
            "validated probabilistic genotyping software (e.g. EuroForMix, STRmix) and "
            "peak-height data"
        ),
        message=(
            "Mixture deconvolution needs a validated probabilistic genotyping model with peak "
            "heights, drop-in/dropout and stutter "
            "parameters; DNAKit does not implement it."
        ),
    )


def infer(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Infer the most probable biological relationship between two profiles from IBD-model
    likelihood ratios.

    input: {'profile_1', 'profile_2', 'frequencies'}. options: relationships (default unrelated,
    parent_child, full_sibling, half_sibling, first_cousin; any of the calculate list), priors
    ({relationship: prior}, default uniform), min_frequency. posterior(r) is proportional to
    prior(r) * LR(r vs unrelated) from the Cotterman model (theta=0, no mutation).
    half_sibling, grandparent_grandchild and uncle_nephew share identical IBD probabilities
    (0.5, 0.5, 0) and are mathematically indistinguishable from autosomal genotypes: their
    likelihood ratios are equal (posteriors differ only through the priors), so whenever two or more
    evaluated relationships are indistinguishable a warning lists them, and most_probable among them
    reflects list order, not evidence.
    """

    return _finish(_M, "infer", input, options, _f.infer)


def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw the per-locus log10 genotype probability / likelihood ratio of a forensic statistic
    as an SVG bar chart.

    input and options as calculate (quantity selects the statistic); path writes the SVG (no
    overwrite unless overwrite=True).
    """

    return _finish(_M, "plot", input, options, _f.plot)


def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Predict phenotype or ancestry from a DNA profile (not implemented: needs validated
    prediction models).
    """

    return candidate_not_implemented(
        _M,
        "predict",
        input,
        options=_get_options(options),
        requires=(
            "validated forensic phenotype/ancestry prediction models and reference data (e.g. "
            "HIrisPlex-S, a reference panel)"
        ),
        message=(
            "Forensic DNA phenotyping and ancestry prediction require validated trained models "
            "and reference populations that are not bundled."
        ),
    )


def read(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Read a genotype profile or an allele-frequency table (path, delimited text, records or
    mappings) into the DNAKit structures.

    Profile columns: locus, allele1, allele2. Frequency columns: locus, allele, frequency.
    options: kind='auto' | 'profile' | 'frequencies'. A mapping with 'profile'/'frequencies'
    keys is normalised too.
    """

    return _finish(_M, "read", input, options, _f.read)


def report(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Generate a Markdown statistical report (per-locus table, combined value, ENFSI verbal
    scale, limitations) for a forensic statistic.

    input and options as calculate (quantity selects the statistic); title; path writes the
    Markdown (no overwrite unless overwrite=True). The statistic itself is included in the
    result.
    """

    return _finish(_M, "report", input, options, _f.report)


def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate forensic profiles and allele-frequency tables (frequency sums, allele counts per
    locus, missing loci/alleles).

    input: mapping of any of 'frequencies', 'profile', 'evidence', 'suspect', 'child', ...
    Errors: frequencies not summing to 1 within tolerance=0.02, more than two alleles at a locus
    (mixture), locus absent from the frequency table. Warnings: single allele listed, allele
    absent from the frequency table. Returns valid=False rather than raising on content
    problems.
    """

    return _finish(_M, "validate", input, options, _f.validate)


__all__ = [
    "calculate",
    "compare",
    "deconvolute",
    "infer",
    "plot",
    "predict",
    "read",
    "report",
    "validate",
]

__candidate_adapter_operations__ = set(
    [
        "calculate",
        "compare",
        "deconvolute",
        "infer",
        "plot",
        "predict",
        "read",
        "report",
        "validate",
    ]
)


from dnakit._verbosity import register_public as _register_public_verbose

_register_public_verbose(globals())
