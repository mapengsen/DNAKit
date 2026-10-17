"""Experimental-design functions: power and sample size, coverage planning, randomised layouts,
assignment, checks and export.

Parameters come from ``input`` (a mapping) and ``options`` (``input`` wins on a name clash for
the numeric functions). Randomisation is reproducible: ``seed`` (default 0) seeds Python's
``random.Random``.  Every function returns the DNAKit candidate envelope.
"""

from __future__ import annotations

from collections.abc import Mapping

from dnakit.gwas._impl_common import finish as _finish

from . import _impl_design as _d

_M = "experimental_design"


def assign(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Randomly assign experimental units to treatment groups by (stratified) permuted-block
    randomisation.

    input: a list of unit ids, an integer number of units, or records (optionally with a
    stratification column). options: groups (names; default G1..Gk from n_groups=2), ratio
    (positive integers per group, default equal), stratify_by (column), seed=0. Within each
    stratum, consecutive blocks of sum(ratio) units receive a random permutation of the group
    pattern.
    """

    return _finish(_M, "assign", input, options, _d.assign)


def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate sample size, power, sequencing coverage (Lander-Waterman) or replicate number
    for common designs.

    Parameters in input mapping or options. quantity = 'sample_size' (default), 'power',
    'coverage' or 'replicates'. test = 'two_sample_t' (default), 'one_sample_t', 'paired_t'
    (effect_size = Cohen's d; exact non-central t), 'two_proportions' (p1, p2;
    pooled-variance normal approximation as R power.prop.test) or 'chi_square' (Cohen's w as
    effect_size, df); alpha=0.05, power=0.8, alternative='two-sided' | 'larger' | 'smaller',
    ratio=1 (n2/n1), n (to compute power). Coverage: genome_size, read_length, n_reads or
    coverage, optional target_fraction, min_depth, overlap_fraction. Replicates: sd, margin,
    confidence=0.95 (smallest n whose t confidence half-width is <= margin).
    """

    return _finish(_M, "calculate", input, options, _d.calculate)


def check(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Check the statistical properties of a design table: replication balance, complete blocks,
    Latin-square property, orthogonality.

    input: records, the full envelope returned by design() (a mapping with 'status' and 'result';
    its 'result' is used, and a status other than 'ok' raises INVALID_DESIGN_ENVELOPE), or the
    inner design() result mapping (with a 'layout' list). The table needs a 'treatment' column
    (options treatment_column), optionally 'block', 'row'/'column'. options: factor_columns (coded numeric columns whose
    pairwise dot products must be 0 for orthogonality).
    """

    return _finish(_M, "check", input, options, _d.check)


def design(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Build a randomised experimental layout: completely randomised, randomised complete block,
    Latin square or factorial.

    input: mapping (merged with options) with design_type = 'crd' (default) | 'rcbd' |
    'latin_square' | 'factorial'; treatments (list) or n_treatments; replicates=2 (CRD
    replicates, or number of blocks for RCBD unless blocks is given); factors ({name: levels})
    for factorial; seed=0. Returns the plot layout and replication counts.
    """

    return _finish(_M, "design", input, options, _d.design)


def export(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Export a design layout as CSV, TSV or JSON text (and optionally a file).

    input: records, the full envelope returned by design() (a mapping with 'status' and 'result';
    its 'result' is used, and a status other than 'ok' raises INVALID_DESIGN_ENVELOPE), or the
    inner design() result mapping (with a 'layout' list). options: format='csv' | 'tsv' | 'json', columns
    (subset/order), path (existing files are not overwritten unless overwrite=True).
    """

    return _finish(_M, "export", input, options, _d.export)


def generate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Generate a standard-order design matrix: full factorial, two-level fractional factorial or
    cyclic Latin square.

    input: mapping (merged with options). kind = 'full_factorial' (factors {name: levels}; first
    factor varies fastest), 'fractional_factorial' (n_factors and generators such as ['D=ABC'];
    coded -1/+1) or 'latin_square' (treatments; cyclic, unrandomised).
    """

    return _finish(_M, "generate", input, options, _d.generate)


def randomize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Randomise the run order of a list of runs (Fisher-Yates shuffle with a fixed seed),
    optionally within groups.

    input: a list of items or records. options: seed=0, within (column; runs are shuffled inside
    each group, groups keep their first-appearance order). Adds run_order (1-based) to every
    run.
    """

    return _finish(_M, "randomize", input, options, _d.randomize)


def simulate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Estimate statistical power by Monte Carlo simulation of a t-test experiment and compare it
    with the exact power.

    input: mapping (merged with options) with n (per group), effect_size (Cohen's d), test =
    'two_sample_t' (default) | 'one_sample_t', alpha=0.05, n_simulations=2000, seed=0. Normal
    data with unit SD, two-sided pooled t-test.
    """

    return _finish(_M, "simulate", input, options, _d.simulate)


def validate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Validate the structure of a design table: required columns present, no empty cells, unique
    plot ids, consistent rows.

    input: records, the full envelope returned by design(), or the inner design() result mapping
    (a non-ok envelope raises INVALID_DESIGN_ENVELOPE). options: required_columns=['treatment'],
    unit_column='plot'. Returns valid=False with the list of errors rather than raising on
    content problems.
    """

    return _finish(_M, "validate", input, options, _d.validate)


__all__ = [
    "assign",
    "calculate",
    "check",
    "design",
    "export",
    "generate",
    "randomize",
    "simulate",
    "validate",
]

__candidate_adapter_operations__ = set(
    [
        "assign",
        "calculate",
        "check",
        "design",
        "export",
        "generate",
        "randomize",
        "simulate",
        "validate",
    ]
)


from dnakit._verbosity import register_public as _register_public_verbose

_register_public_verbose(globals())
