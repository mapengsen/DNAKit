"""scikit-allel wrappers for family-variant rows that no existing dnakit wrapper reaches (group H1).

The earlier ledger called these allel functions variants of ``calculate_fst``,
``calculate_pca``, ``calculate_patterson_d``, ``calculate_garud_h``, the site frequency spectrum
and the VCF reader. None of those wrappers has a parameter that reaches them: they differ in the
output (per variant, per window, with a jackknife error, joint spectra, a field table) or in the
algorithm (randomized PCA, transmission-phasing steps).

Inputs are plain lists of allele counts, genotypes, dosages or haplotypes from local data.
Positions are 1-based and sorted; window bounds are 1-based and inclusive; moving windows count
variants. Statistics that cannot be computed are returned as ``None``. scikit-allel is imported
only when a function runs. Output files are never overwritten.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path
from dnakit.exceptions import ConfigurationError
from dnakit.variants.allel_functions import _alt_counts, _genotypes
from dnakit.variants.allel_review_functions import (
    _MAX_COUNT,
    _allel,
    _bool_array,
    _bounds,
    _clean,
    _counts,
    _fields,
    _fresh_file,
    _haplotypes,
    _index_bound,
    _int_array,
    _numpy,
    _optional_position,
    _optional_window,
    _positions,
    _provenance,
    _quiet,
    _same_rows,
    _sample_size,
    _values,
    _window_rows,
    _window_size,
)

_MAX_ROWS = 100_000
_MAX_BETWEEN = 500
_MAX_COMPONENTS = 50


# ----------------------------------------------------------------------------------------------
# shared input checks
# ----------------------------------------------------------------------------------------------


def _window_options(
    start: object, stop: object, step: object
) -> tuple[int | None, int | None, int | None]:
    left = _optional_position(start, "start")
    right = _optional_position(stop, "stop")
    stride = _optional_window(step, "step")
    _bounds(left, right)
    return left, right, stride


def _moving_options(
    window_size: object, start: object, stop: object, step: object
) -> tuple[int, int, int | None, int | None]:
    size = _window_size(window_size)
    left = _index_bound(start, "start") or 0
    right = _index_bound(stop, "stop")
    stride = _optional_window(step, "step")
    return size, left, right, stride


def _access(is_accessible: object, last_position: int) -> Any:
    if is_accessible is None:
        return None
    access = _bool_array(is_accessible, "is_accessible")
    if access.shape[0] < last_position:
        raise ConfigurationError("is_accessible must cover every position.")
    return access


def _subpopulations(populations: object, samples: int) -> list[list[int]]:
    if isinstance(populations, (str, bytes)) or not isinstance(populations, Sequence):
        raise ConfigurationError("populations must list the sample indices of each population.")
    if len(populations) < 2:
        raise ConfigurationError("populations must list at least two populations.")
    groups: list[list[int]] = []
    seen: set[int] = set()
    for group in populations:
        if isinstance(group, (str, bytes)) or not isinstance(group, Sequence) or not group:
            raise ConfigurationError("each population must be a non-empty list of sample indices.")
        members = []
        for item in group:
            if isinstance(item, bool) or not isinstance(item, int) or not 0 <= item < samples:
                raise ConfigurationError(f"sample indices must be integers in [0, {samples - 1}].")
            members.append(item)
        if seen & set(members):
            raise ConfigurationError("populations must not share samples.")
        seen.update(members)
        groups.append(members)
    return groups


def _four_counts(
    first: object, second: object, third: object, fourth: object
) -> tuple[Any, Any, Any, Any]:
    arrays = (
        _counts(first, "first_counts"),
        _counts(second, "second_counts"),
        _counts(third, "third_counts"),
        _counts(fourth, "fourth_counts"),
    )
    _same_rows(arrays, ["first_counts", "second_counts", "third_counts", "fourth_counts"])
    return arrays


def _count_arrays(allel: Any, *arrays: Any) -> list[Any]:
    return [allel.AlleleCountsArray(array) for array in arrays]


# ----------------------------------------------------------------------------------------------
# windowed and moving statistics
# ----------------------------------------------------------------------------------------------


def calculate_windowed_divergence(
    positions: Sequence[int],
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    is_accessible: Sequence[bool] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Nucleotide divergence (dxy) per base between two populations in windows of base pairs.

    Counts are allele counts per variant for each population. With ``is_accessible`` (one flag per
    base from position 1), only accessible bases count in the denominator. ``calculate_sequence_
    divergence()`` gives one value for a whole region instead.
    """

    allel = _allel()
    pos = _positions(positions)
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    if first.shape[0] != pos.shape[0]:
        raise ConfigurationError("first_counts must have one row per position.")
    size = _window_size(window_size)
    left, right, stride = _window_options(start, stop, step)
    access = _access(is_accessible, int(pos[-1]))
    with _quiet():
        values, windows, bases, sites = allel.windowed_divergence(
            allel.SortedIndex(pos),
            *_count_arrays(allel, first, second),
            size=size,
            start=left,
            stop=right,
            step=stride,
            is_accessible=access,
        )
    rows = _window_rows(windows, "dxy", values, sites)
    for row, count in zip(rows, _numpy().asarray(bases).tolist(), strict=False):
        row["callable_bases"] = int(count)
    return {
        "status": "completed",
        "operation": "calculate_windowed_divergence",
        "window_size": size,
        "step": size if stride is None else stride,
        "windows": rows,
        "provenance": _provenance(),
    }


def calculate_windowed_watterson_theta(
    positions: Sequence[int],
    allele_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    is_accessible: Sequence[bool] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Watterson's theta per base in windows of ``window_size`` base pairs over one contig.

    ``allele_counts`` has one row per variant for one population. With ``is_accessible`` (one flag
    per base from position 1), only accessible bases count in the denominator.
    """

    allel = _allel()
    pos = _positions(positions)
    counts = _counts(allele_counts, "allele_counts")
    if counts.shape[0] != pos.shape[0]:
        raise ConfigurationError("allele_counts must have one row per position.")
    size = _window_size(window_size)
    left, right, stride = _window_options(start, stop, step)
    access = _access(is_accessible, int(pos[-1]))
    with _quiet():
        values, windows, bases, sites = allel.windowed_watterson_theta(
            allel.SortedIndex(pos),
            *_count_arrays(allel, counts),
            size=size,
            start=left,
            stop=right,
            step=stride,
            is_accessible=access,
        )
    rows = _window_rows(windows, "theta_per_base", values, sites)
    for row, count in zip(rows, _numpy().asarray(bases).tolist(), strict=False):
        row["callable_bases"] = int(count)
    return {
        "status": "completed",
        "operation": "calculate_windowed_watterson_theta",
        "window_size": size,
        "step": size if stride is None else stride,
        "windows": rows,
        "provenance": _provenance(),
    }


def calculate_moving_tajimas_d(
    allele_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int = 0,
    stop: int | None = None,
    step: int | None = None,
    min_sites: int = 3,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tajima's D in moving windows of ``window_size`` consecutive variants."""

    allel = _allel()
    counts = _counts(allele_counts, "allele_counts")
    size, left, right, stride = _moving_options(window_size, start, stop, step)
    minimum = _window_size(min_sites, "min_sites")
    with _quiet():
        values = allel.moving_tajima_d(
            allel.AlleleCountsArray(counts),
            size,
            start=left,
            stop=right,
            step=stride,
            min_sites=minimum,
        )
    return {
        "status": "completed",
        "operation": "calculate_moving_tajimas_d",
        "window_size": size,
        "step": size if stride is None else stride,
        "values": _values(values),
        "provenance": _provenance(),
    }


def calculate_weir_cockerham_fst_per_variant(
    genotypes: Sequence[Sequence[Sequence[int]]],
    populations: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Weir-Cockerham variance components and Fst for every variant, using scikit-allel.

    ``genotypes`` has shape (variants, samples, ploidy) with -1 for a missing call, and
    ``populations`` lists the sample indices of each population. For every variant the result
    gives the components ``a`` (among populations), ``b`` (among individuals within populations)
    and ``c`` (within individuals), summed over alleles, and ``fst = a / (a + b + c)``.
    ``calculate_fst()`` returns the genome-wide average only.
    """

    allel = _allel()
    np = _numpy()
    gt = _genotypes(genotypes, "genotypes")
    groups = _subpopulations(populations, int(gt.shape[1]))
    with _quiet():
        a, b, c = allel.weir_cockerham_fst(gt, groups)
    a_sum, b_sum, c_sum = (np.asarray(item, dtype=float).sum(axis=1) for item in (a, b, c))
    with _quiet():
        per_variant = a_sum / (a_sum + b_sum + c_sum)
    total = float(a_sum.sum() + b_sum.sum() + c_sum.sum())
    return {
        "status": "completed",
        "operation": "calculate_weir_cockerham_fst_per_variant",
        "variants": int(gt.shape[0]),
        "populations": [len(group) for group in groups],
        "a": _values(a_sum),
        "b": _values(b_sum),
        "c": _values(c_sum),
        "fst": _values(per_variant),
        "fst_ratio_of_sums": None if total == 0 else _clean(float(a_sum.sum()) / total),
        "provenance": _provenance(),
    }


def calculate_windowed_hudson_fst(
    positions: Sequence[int],
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hudson's Fst (ratio of sums) in windows of ``window_size`` base pairs over one contig."""

    allel = _allel()
    pos = _positions(positions)
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    if first.shape[0] != pos.shape[0]:
        raise ConfigurationError("first_counts must have one row per position.")
    size = _window_size(window_size)
    left, right, stride = _window_options(start, stop, step)
    with _quiet():
        fst, windows, counts = allel.windowed_hudson_fst(
            allel.SortedIndex(pos),
            *_count_arrays(allel, first, second),
            size=size,
            start=left,
            stop=right,
            step=stride,
        )
    return {
        "status": "completed",
        "operation": "calculate_windowed_hudson_fst",
        "window_size": size,
        "step": size if stride is None else stride,
        "windows": _window_rows(windows, "fst", fst, counts),
        "provenance": _provenance(),
    }


def calculate_windowed_weir_cockerham_fst(
    positions: Sequence[int],
    genotypes: Sequence[Sequence[Sequence[int]]],
    populations: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int | None = None,
    stop: int | None = None,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Weir-Cockerham Fst in windows of ``window_size`` base pairs over one contig.

    ``genotypes`` has shape (variants, samples, ploidy) with -1 for a missing call, and
    ``populations`` lists the sample indices of each population.
    """

    allel = _allel()
    pos = _positions(positions)
    gt = _genotypes(genotypes, "genotypes")
    if gt.shape[0] != pos.shape[0]:
        raise ConfigurationError("genotypes must have one row per position.")
    groups = _subpopulations(populations, int(gt.shape[1]))
    size = _window_size(window_size)
    left, right, stride = _window_options(start, stop, step)
    with _quiet():
        fst, windows, counts = allel.windowed_weir_cockerham_fst(
            allel.SortedIndex(pos), gt, groups, size=size, start=left, stop=right, step=stride
        )
    return {
        "status": "completed",
        "operation": "calculate_windowed_weir_cockerham_fst",
        "window_size": size,
        "step": size if stride is None else stride,
        "windows": _window_rows(windows, "fst", fst, counts),
        "provenance": _provenance(),
    }


def calculate_moving_hudson_fst(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int = 0,
    stop: int | None = None,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hudson's Fst in moving windows of ``window_size`` consecutive variants."""

    allel = _allel()
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    size, left, right, stride = _moving_options(window_size, start, stop, step)
    with _quiet():
        values = allel.moving_hudson_fst(
            *_count_arrays(allel, first, second), size=size, start=left, stop=right, step=stride
        )
    return {
        "status": "completed",
        "operation": "calculate_moving_hudson_fst",
        "window_size": size,
        "step": size if stride is None else stride,
        "values": _values(values),
        "provenance": _provenance(),
    }


def calculate_moving_weir_cockerham_fst(
    genotypes: Sequence[Sequence[Sequence[int]]],
    populations: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int = 0,
    stop: int | None = None,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Weir-Cockerham Fst in moving windows of ``window_size`` consecutive variants."""

    allel = _allel()
    gt = _genotypes(genotypes, "genotypes")
    groups = _subpopulations(populations, int(gt.shape[1]))
    size, left, right, stride = _moving_options(window_size, start, stop, step)
    with _quiet():
        values = allel.moving_weir_cockerham_fst(
            gt, groups, size, start=left, stop=right, step=stride
        )
    return {
        "status": "completed",
        "operation": "calculate_moving_weir_cockerham_fst",
        "window_size": size,
        "step": size if stride is None else stride,
        "values": _values(values),
        "provenance": _provenance(),
    }


def calculate_blockwise_patterson_d(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    third_counts: Sequence[Sequence[int]],
    fourth_counts: Sequence[Sequence[int]],
    *,
    block_length: int,
    verbose: bool = True,
) -> dict[str, Any]:
    """Patterson's D (ABBA-BABA) with a block-jackknife standard error and z-score.

    The four populations are A, B, C and D in the order of ``calculate_patterson_d()``; blocks
    hold ``block_length`` consecutive variants. This is also scikit-allel's
    ``average_patterson_d``.
    """

    allel = _allel()
    arrays = _four_counts(first_counts, second_counts, third_counts, fourth_counts)
    blen = _window_size(block_length, "block_length")
    with _quiet():
        d, error, z_score, block_values, jackknife = allel.blockwise_patterson_d(
            *_count_arrays(allel, *arrays), blen
        )
    return {
        "status": "completed",
        "operation": "calculate_blockwise_patterson_d",
        "variants": int(arrays[0].shape[0]),
        "block_length": blen,
        "d": _clean(float(d)),
        "standard_error": _clean(float(error)),
        "z_score": _clean(float(z_score)),
        "block_values": _values(block_values),
        "jackknife_values": _values(jackknife),
        "provenance": _provenance(),
    }


def calculate_moving_patterson_d(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    third_counts: Sequence[Sequence[int]],
    fourth_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int = 0,
    stop: int | None = None,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Patterson's D (ABBA-BABA) in moving windows of ``window_size`` consecutive variants."""

    allel = _allel()
    arrays = _four_counts(first_counts, second_counts, third_counts, fourth_counts)
    size, left, right, stride = _moving_options(window_size, start, stop, step)
    with _quiet():
        values = allel.moving_patterson_d(
            *_count_arrays(allel, *arrays), size=size, start=left, stop=right, step=stride
        )
    return {
        "status": "completed",
        "operation": "calculate_moving_patterson_d",
        "window_size": size,
        "step": size if stride is None else stride,
        "values": _values(values),
        "provenance": _provenance(),
    }


def calculate_moving_garud_h(
    haplotypes: Sequence[Sequence[int]],
    *,
    window_size: int,
    start: int = 0,
    stop: int | None = None,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Garud's H1, H12, H123 and H2/H1 in moving windows of consecutive variants.

    ``haplotypes`` has shape (variants, haplotypes) with alleles 0 and 1. Each window of
    ``window_size`` variants gets its own statistics.
    """

    allel = _allel()
    haps = _haplotypes(haplotypes, "haplotypes")
    size, left, right, stride = _moving_options(window_size, start, stop, step)
    with _quiet():
        h1, h12, h123, h2_h1 = allel.moving_garud_h(haps, size, start=left, stop=right, step=stride)
    return {
        "status": "completed",
        "operation": "calculate_moving_garud_h",
        "window_size": size,
        "step": size if stride is None else stride,
        "h1": _values(h1),
        "h12": _values(h12),
        "h123": _values(h123),
        "h2_h1": _values(h2_h1),
        "provenance": _provenance(),
    }


# ----------------------------------------------------------------------------------------------
# linkage and ordination
# ----------------------------------------------------------------------------------------------


def calculate_rogers_huff_r_between(
    alt_counts_a: Sequence[Sequence[int]],
    alt_counts_b: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rogers-Huff r between every variant of one set and every variant of another.

    Both inputs are alternate-allele dosages (0, 1, 2; -1 for missing) with shape
    (variants, samples) over the same samples. The result has one row per variant of the first
    set and one column per variant of the second; ``calculate_rogers_huff_r()`` pairs the
    variants inside one set instead.
    """

    allel = _allel()
    first = _alt_counts(alt_counts_a, "alt_counts_a")
    second = _alt_counts(alt_counts_b, "alt_counts_b")
    if first.shape[1] != second.shape[1]:
        raise ConfigurationError("alt_counts_a and alt_counts_b must have the same samples.")
    if first.shape[0] > _MAX_BETWEEN or second.shape[0] > _MAX_BETWEEN:
        raise ConfigurationError(f"each set may hold at most {_MAX_BETWEEN} variants.")
    with _quiet():
        matrix = allel.rogers_huff_r_between(first, second)
    return {
        "status": "completed",
        "operation": "calculate_rogers_huff_r_between",
        "variants_a": int(first.shape[0]),
        "variants_b": int(second.shape[0]),
        "r": _values(matrix),
        "provenance": _provenance(),
    }


def calculate_randomized_pca(
    alt_counts: Sequence[Sequence[int]],
    *,
    n_components: int = 10,
    iterated_power: int = 3,
    seed: int = 0,
    scaler: str = "patterson",
    verbose: bool = True,
) -> dict[str, Any]:
    """Randomized principal component analysis of genotype dosages, using scikit-allel.

    ``alt_counts`` has shape (variants, samples) with dosages 0, 1 and 2 and no missing calls.
    Monomorphic variants are dropped. The randomized solver approximates the exact PCA of
    ``calculate_pca()`` and is reproducible through ``seed``. ``scaler`` is ``patterson``
    (default), ``standard`` or ``none``.
    """

    allel = _allel()
    np = _numpy()
    dosage = _alt_counts(alt_counts, "alt_counts")
    if (dosage < 0).any():
        raise ConfigurationError("alt_counts must not contain missing calls (-1).")
    if isinstance(scaler, str) and scaler not in ("patterson", "standard", "none"):
        raise ConfigurationError("scaler must be 'patterson', 'standard' or 'none'.")
    if not isinstance(scaler, str):
        raise ConfigurationError("scaler must be 'patterson', 'standard' or 'none'.")
    components = _window_size(n_components, "n_components")
    if components > _MAX_COMPONENTS:
        raise ConfigurationError(f"n_components may be at most {_MAX_COMPONENTS}.")
    power = _window_size(iterated_power, "iterated_power")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ConfigurationError("seed must be a non-negative integer.")
    total = dosage.sum(axis=1)
    kept = dosage[(total > 0) & (total < 2 * dosage.shape[1])]
    if kept.shape[0] == 0:
        raise ConfigurationError("alt_counts has no polymorphic variants.")
    if components > min(kept.shape):
        raise ConfigurationError(
            "n_components cannot exceed the number of samples or polymorphic variants."
        )
    coordinates, model = allel.randomized_pca(
        kept,
        n_components=components,
        iterated_power=power,
        random_state=seed,
        scaler=None if scaler == "none" else scaler,
    )
    return {
        "status": "completed",
        "operation": "calculate_randomized_pca",
        "n_components": components,
        "seed": seed,
        "scaler": scaler,
        "polymorphic_variants": int(kept.shape[0]),
        "coordinates": [[float(item) for item in row] for row in np.asarray(coordinates)],
        "explained_variance_ratio": [float(item) for item in model.explained_variance_ratio_],
        "provenance": _provenance(),
    }


# ----------------------------------------------------------------------------------------------
# site frequency spectra from allele counts
# ----------------------------------------------------------------------------------------------


def _biallelic_counts(value: object, name: str) -> Any:
    counts = _int_array(value, name, ndim=2, maximum=_MAX_COUNT)
    if counts.shape[1] != 2:
        raise ConfigurationError(f"{name} must have shape (variants, 2): reference and alternate.")
    return counts


def _spectrum_size(counts: Any, given: object, name: str) -> int:
    if given is None:
        size = int(counts.sum(axis=1).max())
        if size < 1:
            raise ConfigurationError(f"{name} has no called alleles.")
        return size
    size = _sample_size(given, name)
    if int(counts.sum(axis=1).max()) > size:
        raise ConfigurationError(f"{name} has a variant with more alleles than the sample size.")
    return size


def _flags(fold: object, scale: object) -> tuple[bool, bool]:
    if not isinstance(fold, bool) or not isinstance(scale, bool):
        raise ConfigurationError("fold and scale must be true or false.")
    return fold, scale


def calculate_site_frequency_spectrum_from_counts(
    allele_counts: Sequence[Sequence[int]],
    *,
    fold: bool = False,
    scale: bool = False,
    sample_size: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Site frequency spectrum from biallelic allele counts, folded and/or scaled.

    ``allele_counts`` has shape (variants, 2): reference and alternate (derived) counts.
    ``spectrum[i]`` counts variants with ``i`` derived alleles of ``sample_size`` chromosomes
    (default: the largest total in the data). ``fold=True`` gives the minor-allele spectrum of
    ``sample_size // 2 + 1`` classes (``sfs_folded``). ``scale=True`` multiplies each class as
    ``scale_site_frequency_spectrum()`` and ``scale_folded_site_frequency_spectrum()`` do
    (``sfs_scaled``, ``sfs_folded_scaled``). The plain form is ``sfs``.
    """

    allel = _allel()
    folded, scaled = _flags(fold, scale)
    counts = _biallelic_counts(allele_counts, "allele_counts")
    size = _spectrum_size(counts, sample_size, "allele_counts")
    if folded:
        spectrum = (
            allel.sfs_folded_scaled(counts, n=size) if scaled else allel.sfs_folded(counts, n=size)
        )
    else:
        derived = counts[:, 1]
        spectrum = allel.sfs_scaled(derived, n=size) if scaled else allel.sfs(derived, n=size)
    return {
        "status": "completed",
        "operation": "calculate_site_frequency_spectrum_from_counts",
        "variants": int(counts.shape[0]),
        "sample_size": size,
        "folded": folded,
        "scaled": scaled,
        "spectrum": _values(spectrum) if scaled else [int(item) for item in spectrum.tolist()],
        "provenance": _provenance(),
    }


def calculate_joint_site_frequency_spectrum_from_counts(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    fold: bool = False,
    scale: bool = False,
    first_sample_size: int | None = None,
    second_sample_size: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Joint site frequency spectrum of two populations from biallelic allele counts.

    Each input has shape (variants, 2): reference and alternate (derived) counts of one
    population, for the same variants. ``spectrum[i][j]`` counts variants with ``i`` derived
    alleles in the first and ``j`` in the second population. ``fold=True`` folds each population
    to its minor allele (``joint_sfs_folded``), which differs from folding the pooled sample in
    ``fold_joint_site_frequency_spectrum()``. ``scale=True`` applies the cell scaling of
    ``scale_joint_site_frequency_spectrum()`` (``joint_sfs_scaled``, ``joint_sfs_folded_scaled``).
    """

    allel = _allel()
    folded, scaled = _flags(fold, scale)
    first = _biallelic_counts(first_counts, "first_counts")
    second = _biallelic_counts(second_counts, "second_counts")
    _same_rows([first, second], ["first_counts", "second_counts"])
    n1 = _spectrum_size(first, first_sample_size, "first_counts")
    n2 = _spectrum_size(second, second_sample_size, "second_counts")
    if folded:
        function = allel.joint_sfs_folded_scaled if scaled else allel.joint_sfs_folded
        spectrum = function(first, second, n1=n1, n2=n2)
    else:
        function = allel.joint_sfs_scaled if scaled else allel.joint_sfs
        spectrum = function(first[:, 1], second[:, 1], n1=n1, n2=n2)
    return {
        "status": "completed",
        "operation": "calculate_joint_site_frequency_spectrum_from_counts",
        "variants": int(first.shape[0]),
        "first_sample_size": n1,
        "second_sample_size": n2,
        "folded": folded,
        "scaled": scaled,
        "spectrum": (
            _values(spectrum)
            if scaled
            else [[int(item) for item in row] for row in spectrum.tolist()]
        ),
        "provenance": _provenance(),
    }


# ----------------------------------------------------------------------------------------------
# transmission phasing steps
# ----------------------------------------------------------------------------------------------


def phase_trio_step_by_transmission(
    genotypes: Sequence[Sequence[Sequence[int]]],
    *,
    step: str = "progeny",
    window_size: int | None = None,
    is_phased: Sequence[Sequence[bool]] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """One step of trio phasing by Mendelian transmission, using scikit-allel.

    ``genotypes`` has shape (variants, 3, 2) for mother, father and child (-1 for missing).
    ``step="progeny"`` phases the child where the parents allow it (``phase_progeny_by_
    transmission``). ``step="parents"`` phases the parents from the phased child within windows
    of ``window_size`` variants (``phase_parents_by_transmission``); it needs the ``genotypes``
    and ``is_phased`` returned by the progeny step. ``phase_trio_by_transmission()`` runs both.
    """

    allel = _allel()
    np = _numpy()
    if step not in ("progeny", "parents"):
        raise ConfigurationError("step must be 'progeny' or 'parents'.")
    gt = _genotypes(genotypes, "genotypes")
    if gt.shape[1] != 3 or gt.shape[2] != 2:
        raise ConfigurationError(
            "genotypes must hold three diploid samples in order: mother, father, child."
        )
    if step == "progeny":
        if window_size is not None or is_phased is not None:
            raise ConfigurationError("window_size and is_phased apply to step='parents' only.")
        phased = allel.phase_progeny_by_transmission(gt)
        size = None
    else:
        if window_size is None or is_phased is None:
            raise ConfigurationError("step='parents' needs window_size and is_phased.")
        size = _window_size(window_size)
        mask = np.asarray(is_phased)
        if mask.dtype != bool or mask.shape != gt.shape[:2]:
            raise ConfigurationError("is_phased must be a boolean array of shape (variants, 3).")
        gt.is_phased = mask
        phased = allel.phase_parents_by_transmission(gt, window_size=size)
    return {
        "status": "completed",
        "operation": "phase_trio_step_by_transmission",
        "step": step,
        "variants": int(gt.shape[0]),
        "window_size": size,
        "genotypes": [
            [[int(allele) for allele in sample] for sample in variant]
            for variant in phased.values.tolist()
        ],
        "is_phased": [[bool(flag) for flag in row] for row in phased.is_phased.tolist()],
        "provenance": _provenance(),
    }


# ----------------------------------------------------------------------------------------------
# VCF field tables
# ----------------------------------------------------------------------------------------------


def _alt_number(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= 20:
        raise ConfigurationError("alt_number must be an integer in [1, 20].")
    return value


def _plain(value: Any) -> Any:
    """A table cell as a JSON value: NaN and missing become None, NumPy scalars plain numbers."""

    np = _numpy()
    if value is None:
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float):
        return None if value != value or value in (float("inf"), float("-inf")) else value
    return value


def read_vcf_field_table(
    path: str,
    *,
    fields: Sequence[str] | None = None,
    alt_number: int = 3,
    max_rows: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read the variant fields of a local VCF into a flat table, using scikit-allel.

    ``fields`` names allel fields such as ``CHROM``, ``POS``, ``ID``, ``REF``, ``ALT``, ``QUAL``,
    ``FILTER_PASS`` or ``variants/DP`` (default: the standard fixed columns). ``ALT`` and
    multi-value INFO fields are spread over ``alt_number`` columns (``ALT_1``, ``ALT_2``, ...).
    Missing values are None. The table has no genotype columns; ``read_vcf()`` returns the same
    records with genotypes. This is the table form of ``vcf_to_dataframe`` and ``vcf_to_recarray``.
    """

    allel = _allel()
    source = input_path(path)
    selected = _fields(fields)
    number = _alt_number(alt_number)
    limit = _window_size(max_rows, "max_rows")
    if limit > _MAX_ROWS:
        raise ConfigurationError(f"max_rows may be at most {_MAX_ROWS}.")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        frame = allel.vcf_to_dataframe(str(source), fields=selected, alt_number=number)
    if len(frame) > limit:
        raise ConfigurationError(f"{source.name} has {len(frame)} records, above max_rows={limit}.")
    columns = [str(name) for name in frame.columns]
    rows = [
        {column: _plain(value) for column, value in zip(columns, record, strict=True)}
        for record in frame.itertuples(index=False, name=None)
    ]
    return {
        "status": "completed",
        "operation": "read_vcf_field_table",
        "source": str(source),
        "alt_number": number,
        "num_rows": len(rows),
        "columns": columns,
        "rows": rows,
        "provenance": _provenance(),
    }


def vcf_to_csv_table(
    source: str,
    output: str,
    *,
    fields: Sequence[str] | None = None,
    alt_number: int = 3,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write the variant fields of a local VCF to a new CSV file, using scikit-allel.

    The columns are those of ``read_vcf_field_table()``: one row per record, with ``ALT`` and
    multi-value INFO fields spread over ``alt_number`` columns. The output must not exist.
    """

    allel = _allel()
    src = input_path(source)
    target = _fresh_file(output, "output", ".csv")
    selected = _fields(fields)
    number = _alt_number(alt_number)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            allel.vcf_to_csv(str(src), str(target), fields=selected, alt_number=number)
    except Exception:
        target.unlink(missing_ok=True)
        raise
    lines = Path(target).read_text(encoding="utf8").splitlines()
    header = lines[0].split(",") if lines else []
    return {
        "status": "completed",
        "operation": "vcf_to_csv_table",
        "source": str(src),
        "output": str(target),
        "alt_number": number,
        "columns": header,
        "num_rows": max(len(lines) - 1, 0),
        "bytes": target.stat().st_size,
        "provenance": _provenance(),
    }


__all__ = [
    "calculate_blockwise_patterson_d",
    "calculate_joint_site_frequency_spectrum_from_counts",
    "calculate_moving_garud_h",
    "calculate_moving_hudson_fst",
    "calculate_moving_patterson_d",
    "calculate_moving_tajimas_d",
    "calculate_moving_weir_cockerham_fst",
    "calculate_randomized_pca",
    "calculate_rogers_huff_r_between",
    "calculate_site_frequency_spectrum_from_counts",
    "calculate_weir_cockerham_fst_per_variant",
    "calculate_windowed_divergence",
    "calculate_windowed_hudson_fst",
    "calculate_windowed_watterson_theta",
    "calculate_windowed_weir_cockerham_fst",
    "phase_trio_step_by_transmission",
    "read_vcf_field_table",
    "vcf_to_csv_table",
]

_register_public_verbose(globals(), tuple(__all__))
