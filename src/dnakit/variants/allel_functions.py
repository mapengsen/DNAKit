"""Offline scikit-allel wrappers for population-genetic statistics on local arrays.

Inputs are plain Python lists (or NumPy arrays) of allele counts, genotype dosages,
phased haplotypes or trio genotypes. Positions are 1-based and must be sorted. Missing
genotype calls are -1. Missing statistic values are returned as ``None``.

scikit-allel is imported only when a function runs; a missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[bio-scikit-allel]"`` hint.
Other allel functions that are not wrapped here are listed in the coverage ledger.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from importlib import metadata
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_ALLEL_EXTRA = "bio-scikit-allel"
_NUMPY_EXTRA = "bio-analysis"
_MAX_VARIANTS = 1_000_000
_MAX_SAMPLES = 100_000
_MAX_WINDOW = 1_000_000_000
_MAX_LD_VARIANTS = 200
_METRICS = ("euclidean", "cityblock", "cosine", "correlation", "hamming")


def _numpy() -> Any:
    return require("numpy", _NUMPY_EXTRA)


def _allel() -> Any:
    return require("allel", _ALLEL_EXTRA)


def _provenance() -> dict[str, str]:
    return {
        "package": f"scikit-allel {metadata.version('scikit-allel')}",
        "reference": "none: offline arrays",
    }


def _array(
    value: object,
    name: str,
    *,
    ndim: int,
    minimum: int = 0,
    maximum: int | None = None,
    allow_missing: bool = False,
    dtype: str = "int64",
) -> Any:
    np = _numpy()
    try:
        arr = np.asarray(value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a rectangular array.") from exc
    if arr.dtype.kind not in "iu":
        raise ConfigurationError(f"{name} must contain integers.")
    if arr.ndim != ndim or arr.size == 0:
        raise ConfigurationError(f"{name} must be a non-empty {ndim}-dimensional array.")
    lowest = -1 if allow_missing else minimum
    if (arr < lowest).any():
        raise ConfigurationError(f"{name} must be at least {lowest}.")
    if maximum is not None and (arr > maximum).any():
        raise ConfigurationError(f"{name} must be at most {maximum}.")
    return arr.astype(dtype)


def _counts(value: object, name: str) -> Any:
    arr = _array(value, name, ndim=2, minimum=0, maximum=_MAX_SAMPLES)
    if arr.shape[1] < 2 or arr.shape[0] > _MAX_VARIANTS:
        raise ConfigurationError(
            f"{name} must have shape (variants, alleles) with at least two alleles."
        )
    return arr


def _positions(value: object, count: int, name: str = "positions") -> Any:
    np = _numpy()
    arr = _array(value, name, ndim=1, minimum=1, maximum=2_000_000_000)
    if arr.shape[0] != count:
        raise ConfigurationError(f"{name} must have one entry per variant ({count}).")
    if (np.diff(arr) < 0).any():
        raise ConfigurationError(f"{name} must be sorted in ascending order.")
    return arr


def _same_variants(first: Any, second: Any, name: str) -> None:
    if first.shape[0] != second.shape[0]:
        raise ConfigurationError(f"{name} must describe the same number of variants.")


def _genotypes(value: object, name: str) -> Any:
    allel = _allel()
    arr = _array(value, name, ndim=3, minimum=0, maximum=127, allow_missing=True, dtype="int8")
    if arr.shape[2] < 1 or arr.shape[1] > _MAX_SAMPLES:
        raise ConfigurationError(f"{name} must have shape (variants, samples, ploidy).")
    return allel.GenotypeArray(arr)


def _haplotypes(value: object, name: str) -> Any:
    allel = _allel()
    arr = _array(value, name, ndim=2, minimum=0, maximum=1, dtype="int8")
    if arr.shape[1] < 2:
        raise ConfigurationError(f"{name} must have at least two haplotypes.")
    return allel.HaplotypeArray(arr)


def _alt_counts(value: object, name: str) -> Any:
    """Alternate-allele dosages with shape (variants, samples); -1 marks a missing call."""
    arr = _array(value, name, ndim=2, minimum=0, maximum=2, allow_missing=True, dtype="int8")
    if arr.shape[1] < 2:
        raise ConfigurationError(f"{name} must have at least two samples.")
    return arr


def _clean(item: Any) -> Any:
    if isinstance(item, list):
        return [_clean(child) for child in item]
    if item is None:
        return None
    if isinstance(item, float):
        return None if math.isnan(item) or math.isinf(item) else round(item, 10)
    return item


def _listed(value: Any) -> Any:
    np = _numpy()
    return _clean(np.asarray(value, dtype=float).tolist())


def _ints(value: Any) -> Any:
    np = _numpy()
    return np.asarray(value).astype(int).tolist()


def _window_rows(values: Any, windows: Any, counts: Any, label: str) -> list[dict[str, Any]]:
    np = _numpy()
    rows = []
    for (start, stop), site_count, value in zip(
        np.asarray(windows).tolist(),
        np.asarray(counts).tolist(),
        np.asarray(values, dtype=float).tolist(),
        strict=False,
    ):
        rows.append(
            {
                "start": int(start),
                "stop": int(stop),
                "sites": int(site_count),
                label: _clean(float(value)),
            }
        )
    return rows


def _window_size(value: object, name: str = "window_size") -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= _MAX_WINDOW:
        raise ConfigurationError(f"{name} must be a positive integer.")
    return value


def _optional_step(value: object, name: str = "step") -> int | None:
    if value is None:
        return None
    return _window_size(value, name)


def _bounded_float(value: object, name: str, low: float, high: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not low <= value <= high:
        raise ConfigurationError(f"{name} must be a number in [{low}, {high}].")
    return float(value)


def calculate_windowed_tajimas_d(
    positions: Sequence[int],
    allele_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    step: int | None = None,
    min_sites: int = 3,
    verbose: bool = True,
) -> dict[str, Any]:
    """Tajima's D in windows of ``window_size`` base pairs over one contig."""

    allel = _allel()
    counts = _counts(allele_counts, "allele_counts")
    pos = _positions(positions, counts.shape[0])
    size = _window_size(window_size)
    stride = _optional_step(step)
    minimum = _window_size(min_sites, "min_sites")
    values, windows, sites = allel.windowed_tajima_d(
        pos, allel.AlleleCountsArray(counts), size=size, step=stride, min_sites=minimum
    )
    return {
        "status": "completed",
        "operation": "calculate_windowed_tajimas_d",
        "window_size": size,
        "step": size if stride is None else stride,
        "windows": _window_rows(values, windows, sites, "tajima_d"),
        "provenance": _provenance(),
    }


def calculate_windowed_nucleotide_diversity(
    positions: Sequence[int],
    allele_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    step: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Nucleotide diversity (pi) in windows of ``window_size`` base pairs over one contig."""

    allel = _allel()
    counts = _counts(allele_counts, "allele_counts")
    pos = _positions(positions, counts.shape[0])
    size = _window_size(window_size)
    stride = _optional_step(step)
    diversity, windows, bases, sites = allel.windowed_diversity(
        pos, allel.AlleleCountsArray(counts), size=size, step=stride
    )
    rows = []
    for (start, stop), n_bases, n_sites, value in zip(
        windows.tolist(), bases.tolist(), sites.tolist(), _listed(diversity), strict=False
    ):
        rows.append(
            {
                "start": int(start),
                "stop": int(stop),
                "callable_bases": int(n_bases),
                "sites": int(n_sites),
                "pi": value,
            }
        )
    return {
        "status": "completed",
        "operation": "calculate_windowed_nucleotide_diversity",
        "window_size": size,
        "step": size if stride is None else stride,
        "windows": rows,
        "provenance": _provenance(),
    }


def calculate_mean_pairwise_difference(
    allele_counts: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mean number of pairwise differences per variant, from allele counts."""

    allel = _allel()
    counts = _counts(allele_counts, "allele_counts")
    values = allel.mean_pairwise_difference(allel.AlleleCountsArray(counts))
    return {
        "status": "completed",
        "operation": "calculate_mean_pairwise_difference",
        "variants": counts.shape[0],
        "values": _listed(values),
        "provenance": _provenance(),
    }


def calculate_hudson_fst_per_site(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Hudson's Fst numerator and denominator per variant, with the ratio of their sums."""

    allel = _allel()
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_variants(first, second, "second_counts")
    numerator, denominator = allel.hudson_fst(
        allel.AlleleCountsArray(first), allel.AlleleCountsArray(second)
    )
    total_den = float(denominator.sum())
    return {
        "status": "completed",
        "operation": "calculate_hudson_fst_per_site",
        "variants": first.shape[0],
        "numerator": _listed(numerator),
        "denominator": _listed(denominator),
        "fst": None if total_den == 0 else _clean(float(numerator.sum()) / total_den),
        "provenance": _provenance(),
    }


def calculate_patterson_d(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    third_counts: Sequence[Sequence[int]],
    fourth_counts: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Patterson's D statistic D(A, B; C, D) (ABBA-BABA) from four population allele counts."""

    allel = _allel()
    arrays = [
        _counts(first_counts, "first_counts"),
        _counts(second_counts, "second_counts"),
        _counts(third_counts, "third_counts"),
        _counts(fourth_counts, "fourth_counts"),
    ]
    for label, array in zip(
        ("second_counts", "third_counts", "fourth_counts"), arrays[1:], strict=False
    ):
        _same_variants(arrays[0], array, label)
    wrapped = [allel.AlleleCountsArray(array) for array in arrays]
    numerator, denominator = allel.patterson_d(*wrapped)
    total_den = float(denominator.sum())
    return {
        "status": "completed",
        "operation": "calculate_patterson_d",
        "variants": arrays[0].shape[0],
        "numerator": _listed(numerator),
        "denominator": _listed(denominator),
        "d": None if total_den == 0 else _clean(float(numerator.sum()) / total_den),
        "provenance": _provenance(),
    }


def calculate_expected_heterozygosity(
    allele_frequencies: Sequence[Sequence[float]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Expected heterozygosity 1 - sum(p^2) per variant, from allele frequencies in each row."""

    allel = _allel()
    np = _numpy()
    try:
        freqs = np.asarray(allele_frequencies, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError("allele_frequencies must be a rectangular numeric array.") from exc
    if freqs.ndim != 2 or freqs.size == 0 or freqs.shape[1] < 2:
        raise ConfigurationError(
            "allele_frequencies must have shape (variants, alleles) with at least two alleles."
        )
    if not np.isfinite(freqs).all() or (freqs < 0).any() or (freqs > 1).any():
        raise ConfigurationError("allele_frequencies must be finite values in [0, 1].")
    # allel's heterozygosity_expected uses its ploidy argument as an exponent, so only the
    # diploid form 1 - sum(p^2) is exposed here.
    values = allel.heterozygosity_expected(freqs, 2)
    return {
        "status": "completed",
        "operation": "calculate_expected_heterozygosity",
        "ploidy": 2,
        "values": _listed(values),
        "provenance": _provenance(),
    }


def calculate_inbreeding_coefficient(
    genotypes: Sequence[Sequence[Sequence[int]]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Inbreeding coefficient F = 1 - Ho / He for each variant, from diploid genotypes."""

    allel = _allel()
    gt = _genotypes(genotypes, "genotypes")
    values = allel.inbreeding_coefficient(gt)
    return {
        "status": "completed",
        "operation": "calculate_inbreeding_coefficient",
        "variants": int(gt.shape[0]),
        "values": _listed(values),
        "provenance": _provenance(),
    }


def calculate_rogers_huff_r(
    alt_counts: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Rogers-Huff r (Pearson correlation of genotype dosages) for every pair of variants."""

    allel = _allel()
    dosage = _alt_counts(alt_counts, "alt_counts")
    variants = int(dosage.shape[0])
    if variants < 2:
        raise ConfigurationError("alt_counts must contain at least two variants.")
    if variants > _MAX_LD_VARIANTS:
        raise ConfigurationError(
            f"alt_counts may hold at most {_MAX_LD_VARIANTS} variants for an all-pairs table."
        )
    condensed = allel.rogers_huff_r(dosage)
    pairs = []
    index = 0
    for first in range(variants):
        for second in range(first + 1, variants):
            pairs.append(
                {"variant_a": first, "variant_b": second, "r": _clean(float(condensed[index]))}
            )
            index += 1
    return {
        "status": "completed",
        "operation": "calculate_rogers_huff_r",
        "variants": variants,
        "pairs": pairs,
        "provenance": _provenance(),
    }


def select_unlinked_variants(
    alt_counts: Sequence[Sequence[int]],
    *,
    size: int = 100,
    step: int = 20,
    threshold: float = 0.1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Keep variants in approximate linkage equilibrium by sliding-window r-squared pruning."""

    allel = _allel()
    dosage = _alt_counts(alt_counts, "alt_counts")
    window = _window_size(size, "size")
    stride = _window_size(step, "step")
    limit = _bounded_float(threshold, "threshold", 0.0, 1.0)
    mask = allel.locate_unlinked(dosage, size=window, step=stride, threshold=limit)
    kept = [int(index) for index, flag in enumerate(mask.tolist()) if flag]
    return {
        "status": "completed",
        "operation": "select_unlinked_variants",
        "variants": int(dosage.shape[0]),
        "size": window,
        "step": stride,
        "threshold": limit,
        "kept_indices": kept,
        "kept_count": len(kept),
        "provenance": _provenance(),
    }


def calculate_pairwise_distance(
    alt_counts: Sequence[Sequence[int]],
    *,
    metric: str = "euclidean",
    verbose: bool = True,
) -> dict[str, Any]:
    """Square matrix of distances between samples, from alternate-allele dosages.

    Input rows are variants and columns are samples.
    """

    allel = _allel()
    np = _numpy()
    dosage = _alt_counts(alt_counts, "alt_counts")
    if metric not in _METRICS:
        raise ConfigurationError(f"metric must be one of: {', '.join(_METRICS)}.")
    if dosage.shape[0] < 2:
        raise ConfigurationError("alt_counts must contain at least two variants.")
    if (dosage < 0).any():
        raise ConfigurationError(
            "alt_counts cannot contain missing calls for distance computation."
        )
    condensed = np.asarray(allel.pairwise_distance(dosage.astype(float), metric))
    samples = int(dosage.shape[1])
    square = np.zeros((samples, samples))
    upper = np.triu_indices(samples, k=1)
    square[upper] = condensed
    square = square + square.T
    return {
        "status": "completed",
        "operation": "calculate_pairwise_distance",
        "metric": metric,
        "samples": samples,
        "matrix": _listed(square),
        "provenance": _provenance(),
    }


def calculate_ihs(
    haplotypes: Sequence[Sequence[int]],
    positions: Sequence[int],
    *,
    min_ehh: float = 0.05,
    min_maf: float = 0.05,
    gap_scale: int = 20_000,
    max_gap: int = 200_000,
    include_edges: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Unstandardized integrated haplotype score (iHS) per variant from phased haplotypes."""

    allel = _allel()
    haps = _haplotypes(haplotypes, "haplotypes")
    pos = _positions(positions, haps.shape[0])
    values = allel.ihs(
        haps,
        pos,
        min_ehh=_bounded_float(min_ehh, "min_ehh", 0.0, 1.0),
        min_maf=_bounded_float(min_maf, "min_maf", 0.0, 1.0),
        gap_scale=_window_size(gap_scale, "gap_scale"),
        max_gap=_window_size(max_gap, "max_gap"),
        include_edges=bool(include_edges),
        use_threads=False,
    )
    return {
        "status": "completed",
        "operation": "calculate_ihs",
        "variants": int(haps.shape[0]),
        "ihs": _listed(values),
        "provenance": _provenance(),
    }


def calculate_xpehh(
    haplotypes_a: Sequence[Sequence[int]],
    haplotypes_b: Sequence[Sequence[int]],
    positions: Sequence[int],
    *,
    min_ehh: float = 0.05,
    gap_scale: int = 20_000,
    max_gap: int = 200_000,
    include_edges: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Unstandardized cross-population EHH score (XP-EHH) per variant for two haplotype sets."""

    allel = _allel()
    first = _haplotypes(haplotypes_a, "haplotypes_a")
    second = _haplotypes(haplotypes_b, "haplotypes_b")
    _same_variants(first, second, "haplotypes_b")
    pos = _positions(positions, first.shape[0])
    values = allel.xpehh(
        first,
        second,
        pos,
        min_ehh=_bounded_float(min_ehh, "min_ehh", 0.0, 1.0),
        gap_scale=_window_size(gap_scale, "gap_scale"),
        max_gap=_window_size(max_gap, "max_gap"),
        include_edges=bool(include_edges),
        use_threads=False,
    )
    return {
        "status": "completed",
        "operation": "calculate_xpehh",
        "variants": int(first.shape[0]),
        "xpehh": _listed(values),
        "provenance": _provenance(),
    }


def calculate_nsl(
    haplotypes: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Unstandardized number of segregating sites by length (nSL) per variant."""

    allel = _allel()
    haps = _haplotypes(haplotypes, "haplotypes")
    values = allel.nsl(haps, use_threads=False)
    return {
        "status": "completed",
        "operation": "calculate_nsl",
        "variants": int(haps.shape[0]),
        "nsl": _listed(values),
        "provenance": _provenance(),
    }


def calculate_garud_h(
    haplotypes: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Garud's H1, H12, H123 and H2/H1 statistics for soft sweeps from phased haplotypes."""

    allel = _allel()
    haps = _haplotypes(haplotypes, "haplotypes")
    h1, h12, h123, h2_h1 = (float(value) for value in allel.garud_h(haps))
    return {
        "status": "completed",
        "operation": "calculate_garud_h",
        "haplotypes": int(haps.shape[1]),
        "h1": _clean(h1),
        "h12": _clean(h12),
        "h123": _clean(h123),
        "h2_h1": _clean(h2_h1),
        "provenance": _provenance(),
    }


def calculate_ehh_decay(
    haplotypes: Sequence[Sequence[int]],
    *,
    truncate: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Decay of extended haplotype homozygosity (EHH) away from the first variant."""

    allel = _allel()
    haps = _haplotypes(haplotypes, "haplotypes")
    values = allel.ehh_decay(haps, truncate=bool(truncate))
    return {
        "status": "completed",
        "operation": "calculate_ehh_decay",
        "variants": int(haps.shape[0]),
        "ehh": _listed(values),
        "provenance": _provenance(),
    }


def find_mendel_errors(
    parents: Sequence[Sequence[Sequence[int]]],
    progeny: Sequence[Sequence[Sequence[int]]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count Mendelian inconsistencies per variant for each progeny of a mother-father pair."""

    allel = _allel()
    parent_gt = _genotypes(parents, "parents")
    progeny_gt = _genotypes(progeny, "progeny")
    if parent_gt.shape[1] != 2:
        raise ConfigurationError("parents must hold exactly two samples: mother and father.")
    if parent_gt.shape[0] != progeny_gt.shape[0]:
        raise ConfigurationError("parents and progeny must describe the same variants.")
    if parent_gt.shape[2] != 2 or progeny_gt.shape[2] != 2:
        raise ConfigurationError("Mendelian checks need diploid genotypes.")
    errors = allel.mendel_errors(parent_gt, progeny_gt)
    per_variant = [[int(count) for count in row] for row in errors.tolist()]
    return {
        "status": "completed",
        "operation": "find_mendel_errors",
        "variants": int(parent_gt.shape[0]),
        "progeny": int(progeny_gt.shape[1]),
        "errors": per_variant,
        "total_errors": int(sum(sum(row) for row in per_variant)),
        "provenance": _provenance(),
    }


def phase_trio_by_transmission(
    genotypes: Sequence[Sequence[Sequence[int]]],
    window_size: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Phase a mother-father-child trio by Mendelian transmission within windows of variants."""

    allel = _allel()
    gt = _genotypes(genotypes, "genotypes")
    if gt.shape[1] != 3:
        raise ConfigurationError(
            "genotypes must hold three samples in order: mother, father, child."
        )
    if gt.shape[2] != 2:
        raise ConfigurationError("Transmission phasing needs diploid genotypes.")
    size = _window_size(window_size)
    phased = allel.phase_by_transmission(gt, window_size=size)
    values = [
        [[int(allele) for allele in sample] for sample in variant]
        for variant in phased.values.tolist()
    ]
    flags = [[bool(flag) for flag in row] for row in phased.is_phased.tolist()]
    return {
        "status": "completed",
        "operation": "phase_trio_by_transmission",
        "variants": int(gt.shape[0]),
        "window_size": size,
        "genotypes": values,
        "phased": flags,
        "provenance": _provenance(),
    }


def call_roh_hmm(
    genotypes: Sequence[Sequence[int]],
    positions: Sequence[int],
    *,
    phet_roh: float = 0.001,
    transition: float = 1e-6,
    min_roh: int = 0,
    contig_size: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Runs of homozygosity for one diploid individual from its genotypes (variants x ploidy)."""

    allel = _allel()
    np = _numpy()
    calls = _array(
        genotypes, "genotypes", ndim=2, minimum=0, maximum=1, allow_missing=True, dtype="int8"
    )
    if calls.shape[1] != 2:
        raise ConfigurationError("genotypes must have two alleles per variant (diploid).")
    pos = _positions(positions, calls.shape[0])
    het = _bounded_float(phet_roh, "phet_roh", 0.0, 1.0)
    jump = _bounded_float(transition, "transition", 0.0, 1.0)
    if isinstance(min_roh, bool) or not isinstance(min_roh, int) or min_roh < 0:
        raise ConfigurationError("min_roh must be a non-negative integer.")
    # allel indexes its observation vector by position, so the contig must reach the last variant.
    size = int(pos[-1]) if contig_size is None else _window_size(contig_size, "contig_size")
    if size < int(pos[-1]):
        raise ConfigurationError("contig_size must be at least the last position.")
    runs, froh = allel.roh_mhmm(
        allel.GenotypeVector(calls),
        pos,
        phet_roh=het,
        transition=jump,
        min_roh=min_roh,
        contig_size=size,
    )
    rows = []
    for record in runs.to_dict("records"):
        rows.append(
            {
                "start": int(record["start"]),
                "stop": int(record["stop"]),
                "length": int(record["length"]),
                "is_marginal": bool(record["is_marginal"]),
            }
        )
    rows.sort(key=lambda row: (row["start"], row["stop"]))
    return {
        "status": "completed",
        "operation": "call_roh_hmm",
        "variants": int(calls.shape[0]),
        "runs": rows,
        "froh": _clean(float(np.asarray(froh))),
        "provenance": _provenance(),
    }


def calculate_sequence_divergence(
    positions: Sequence[int],
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    *,
    start: int | None = None,
    stop: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Nucleotide divergence (dxy) between two populations over a region of the contig."""

    allel = _allel()
    first = _counts(first_counts, "first_counts")
    second = _counts(second_counts, "second_counts")
    _same_variants(first, second, "second_counts")
    pos = _positions(positions, first.shape[0])
    left = None if start is None else _window_size(start, "start")
    right = None if stop is None else _window_size(stop, "stop")
    if left is not None and right is not None and right < left:
        raise ConfigurationError("stop must be at least start.")
    value = allel.sequence_divergence(
        pos,
        allel.AlleleCountsArray(first),
        allel.AlleleCountsArray(second),
        start=left,
        stop=right,
    )
    return {
        "status": "completed",
        "operation": "calculate_sequence_divergence",
        "variants": int(first.shape[0]),
        "start": int(pos[0]) if left is None else left,
        "stop": int(pos[-1]) if right is None else right,
        "dxy": _clean(float(value)),
        "provenance": _provenance(),
    }


def calculate_population_branching_statistic(
    first_counts: Sequence[Sequence[int]],
    second_counts: Sequence[Sequence[int]],
    third_counts: Sequence[Sequence[int]],
    *,
    window_size: int,
    window_step: int | None = None,
    normed: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Population branching statistic (PBS) in windows of ``window_size`` variants."""

    allel = _allel()
    arrays = [
        _counts(first_counts, "first_counts"),
        _counts(second_counts, "second_counts"),
        _counts(third_counts, "third_counts"),
    ]
    for label, array in zip(("second_counts", "third_counts"), arrays[1:], strict=False):
        _same_variants(arrays[0], array, label)
    size = _window_size(window_size)
    stride = _optional_step(window_step, "window_step")
    values = allel.pbs(
        *[allel.AlleleCountsArray(array) for array in arrays],
        window_size=size,
        window_step=stride,
        normed=bool(normed),
    )
    return {
        "status": "completed",
        "operation": "calculate_population_branching_statistic",
        "variants": int(arrays[0].shape[0]),
        "window_size": size,
        "window_step": size if stride is None else stride,
        "normed": bool(normed),
        "values": _listed(values),
        "provenance": _provenance(),
    }


__all__ = [
    "calculate_ehh_decay",
    "calculate_expected_heterozygosity",
    "calculate_garud_h",
    "calculate_hudson_fst_per_site",
    "calculate_ihs",
    "calculate_inbreeding_coefficient",
    "calculate_mean_pairwise_difference",
    "calculate_nsl",
    "calculate_pairwise_distance",
    "calculate_patterson_d",
    "calculate_population_branching_statistic",
    "calculate_rogers_huff_r",
    "calculate_sequence_divergence",
    "calculate_windowed_nucleotide_diversity",
    "calculate_windowed_tajimas_d",
    "calculate_xpehh",
    "call_roh_hmm",
    "find_mendel_errors",
    "phase_trio_by_transmission",
    "select_unlinked_variants",
]

_register_public_verbose(globals(), tuple(__all__))
