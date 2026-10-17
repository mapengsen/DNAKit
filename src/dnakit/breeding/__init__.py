"""Breeding functions: breeder's equation, pedigree/genomic relationships, GBLUP/ridge prediction,
selection, optimum-contribution selection and forward breeding simulation (native,
deterministic).
"""

from __future__ import annotations

from collections.abc import Mapping

from dnakit.gwas._impl_common import finish as _finish

from . import _impl_breeding as _b

_M = "breeding"


def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate a quantitative-genetics quantity: breeder's equation response, heritability,
    pedigree or genomic relationships.

    options: quantity = 'breeders_equation' (default; input mapping with heritability and
    selection_differential, or selection_intensity / proportion_selected with phenotypic_sd,
    optional generation_interval: R = h2*S, S = i*sigma_p, truncation intensity i = phi(z)/p),
    'heritability' (input: parent+offspring lists regressed (relationship='midparent' or
    'single_parent'), or response+selection_differential, or
    genetic_variance+phenotypic_variance), 'pedigree_relationship' (input: records id/sire/dam;
    Henderson tabular A matrix and inbreeding F = A_ii - 1) or 'genomic_relationship' (input:
    genotype matrix; VanRaden GRM, numpy needed).
    """

    return _finish(_M, "calculate", input, options, _b.calculate)


def fit(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Fit a genomic ridge-regression (SNP-BLUP, equivalent to GBLUP) model of phenotype on
    marker dosages.

    input: {'genotypes': samples x markers dosages, 'phenotypes': [...]} (None phenotypes are
    dropped). options: ridge_lambda (penalty on centred dosages) or h2 (then lambda =
    2*sum(p*q)*(1-h2)/h2 with p from the training sample). Returns marker effects,
    intercepts, training predictions and breeding values (numpy needed; not_implemented
    without it).
    """

    return _finish(_M, "fit", input, options, _b.fit)


def optimize(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Optimise parental contributions (optimum-contribution selection): maximise mean breeding
    value minus a coancestry penalty.

    input: {'candidates': records with id and ebv, 'relationship': {'ids': [...], 'matrix':
    [[...]]}} or {'candidates': ..., 'pedigree': records id/sire/dam}. Solves max c'u -
    (penalty/2) c'Ac with sum(c)=1, c>=0 by an active-set method. options: penalty=1.0,
    trait='ebv'. Returns contributions, mean EBV and group coancestry (c'Ac/2).
    """

    return _finish(_M, "optimize", input, options, _b.optimize)


def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Draw breeding progress (values per generation) as an SVG line chart.

    input: {'series': {name: [values per generation]}} or a breeding.simulate result
    (options['metric']='mean_genetic_value' selects the plotted column). options: path (write
    the SVG; no overwrite unless overwrite=True).
    """

    return _finish(_M, "plot", input, options, _b.plot)


def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Predict genomic estimated breeding values of new individuals with a fitted ridge/GBLUP
    model.

    input: {'model': <result of breeding.fit>, 'genotypes': new dosages} or {'train_genotypes',
    'phenotypes', 'genotypes'} (fits first; options ridge_lambda or h2 as in fit; numpy needed).
    Missing dosages are replaced by the training allele frequency (model path) or training mean
    (refit path).
    """

    return _finish(_M, "predict", input, options, _b.predict)


def select(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Select breeding candidates by truncation on a trait/EBV or by a Smith-Hazel selection
    index.

    input: records with id and trait columns. options: method='truncation' (default;
    trait='ebv') or 'index' (traits, economic_weights, genetic_covariance G,
    phenotypic_covariance P; weights b = P^-1 G a); exactly one of n_select or proportion;
    higher_is_better=True. Returns the selected ids, ranking and the selection differential.
    """

    return _finish(_M, "select", input, options, _b.select)


def simulate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Simulate recurrent truncation selection on a polygenic trait (unlinked additive loci) and
    report response per generation.

    input: None or a mapping of the options (overrides options). options: n_individuals=1000,
    n_loci=100, heritability=0.5 (base population), selection_proportion=0.2, generations=5,
    seed=0. Random mating of the selected parents with free recombination; reports mean genetic
    value, variances, selection differential, observed and breeder's-equation-predicted response
    per generation.
    """

    return _finish(_M, "simulate", input, options, _b.simulate)


__all__ = ["calculate", "fit", "optimize", "plot", "predict", "select", "simulate"]

__candidate_adapter_operations__ = set(
    ["calculate", "fit", "optimize", "plot", "predict", "select", "simulate"]
)


from dnakit._verbosity import register_public as _register_public_verbose

_register_public_verbose(globals())
