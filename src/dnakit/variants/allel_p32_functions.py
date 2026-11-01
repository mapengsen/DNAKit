"""scikit-allel genotype, allele-count, position-index and state-run operations (integration pack P32).

Inputs are plain Python lists (or NumPy arrays) from local data: genotype arrays are
``variants x samples x ploidy`` integer lists with ``-1`` for a missing call, haplotype arrays are
``variants x haplotypes``, allele-count arrays are ``variants x alleles``. Positions are 1-based and
sorted. Results are JSON-serialisable dictionaries. One function covers all the scikit-allel methods
that only differ by the kind of call, the allele or the array type, and the docstring names the values
of the selecting parameter.

scikit-allel is imported when a function runs; a missing package raises ``BackendUnavailableError``
with the ``pip install "dnakit[bio-scikit-allel]"`` hint. Tested with scikit-allel 1.3.13.
"""

from __future__ import annotations

import os
from collections.abc import Mapping, Sequence
from importlib import metadata
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError

_ALLEL_EXTRA = "bio-scikit-allel"
_NUMPY_EXTRA = "bio-analysis"
_MAX_ITEMS = 5_000_000
_MAX_ALLELE = 1000
_MAX_POSITION = 2_000_000_000
_MAX_COUNT = 2_000_000_000


def _numpy() -> Any:
    return require("numpy", _NUMPY_EXTRA)


def _allel() -> Any:
    return require("allel", _ALLEL_EXTRA)


def _provenance() -> dict[str, str]:
    return {
        "package": f"scikit-allel {metadata.version('scikit-allel')}",
        "reference": "none: offline arrays or local files",
    }


def _choice(value: object, name: str, allowed: Sequence[str]) -> str:
    if not isinstance(value, str) or value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}.")
    return value


def _int_option(value: object, name: str, *, minimum: int = 0, maximum: int = _MAX_ALLELE) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _int_array(value: object, name: str, *, ndim: int, minimum: int = -1, maximum: int = _MAX_ALLELE) -> Any:
    np = _numpy()
    try:
        array = np.asarray(value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{name} must be a rectangular array.") from exc
    if array.dtype.kind not in "iu" or array.ndim != ndim or array.size == 0:
        raise ConfigurationError(f"{name} must be a non-empty {ndim}-dimensional integer array.")
    if array.size > _MAX_ITEMS:
        raise ConfigurationError(f"{name} exceeds {_MAX_ITEMS} values.")
    if (array < minimum).any() or (array > maximum).any():
        raise ConfigurationError(f"{name} must hold integers in [{minimum}, {maximum}].")
    return array.astype(np.int64)


def _genotypes(value: object, name: str = "genotypes") -> Any:
    return _allel().GenotypeArray(_int_array(value, name, ndim=3))


def _haplotypes(value: object, name: str = "haplotypes") -> Any:
    return _allel().HaplotypeArray(_int_array(value, name, ndim=2))


def _allele_counts(value: object, name: str = "allele_counts") -> Any:
    return _allel().AlleleCountsArray(_int_array(value, name, ndim=2, minimum=0, maximum=_MAX_COUNT))


def _positions(value: object, name: str) -> Any:
    np = _numpy()
    array = _int_array(value, name, ndim=1, minimum=1, maximum=_MAX_POSITION)
    if (np.diff(array) < 0).any():
        raise ConfigurationError(f"{name} must be sorted in ascending order.")
    return array


def _plain(value: Any) -> Any:
    np = _numpy()
    if isinstance(value, np.ndarray):
        if value.dtype.kind == "f" and not np.isfinite(value).all():
            return np.where(np.isfinite(value), value, None).tolist()  # NaN and infinity become null
        if value.dtype.kind in "SU":
            return np.char.decode(value.astype("S"), "ascii").tolist() if value.dtype.kind == "S" else value.tolist()
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def _result(operation: str, **fields: Any) -> dict[str, Any]:
    return {"status": "completed", "operation": operation, **fields, "provenance": _provenance()}


def _new_file(value: object, name: str) -> Path:
    if not isinstance(value, (str, os.PathLike)) or not str(value):
        raise ConfigurationError(f"{name} must be a path.")
    target = Path(value).expanduser().resolve()
    if target.exists():
        raise ConfigurationError(f"{name} must not exist yet: {target}")
    if not target.parent.is_dir():
        raise ConfigurationError(f"The folder for {name} does not exist: {target.parent}")
    return target


def write_sites_vcf(
    output_vcf: str,
    variants: Sequence[Mapping[str, Any]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a sites-only VCF (no sample columns) from variant records with ``allel.write_vcf``.

    Each record has ``chrom``, 1-based ``pos``, ``ref`` and ``alt`` (a string or a list of strings) and
    optionally ``id``, ``qual`` (missing is written as ``.``) and ``filters`` (a list of failed filter names; none means PASS). Covers ``write_vcf`` for a variant table
    (``VariantTable.to_vcf`` and ``VariantChunkedTable.to_vcf`` build the same callset). An existing file is refused.
    """

    np = _numpy()
    allel = _allel()
    target = _new_file(output_vcf, "output_vcf")
    if isinstance(variants, (str, bytes)) or not isinstance(variants, Sequence) or not variants:
        raise ConfigurationError("variants must be a non-empty list of records.")
    chrom, pos, ref, ids, qual, failed, alts = [], [], [], [], [], [], []
    for record in variants:
        if not isinstance(record, Mapping) or not {"chrom", "pos", "ref", "alt"} <= set(record):
            raise ConfigurationError("Each variant needs chrom, pos, ref and alt.")
        alt = record["alt"]
        alt_list = [alt] if isinstance(alt, str) else list(alt)
        if not alt_list or not all(isinstance(item, str) and item for item in alt_list):
            raise ConfigurationError("alt must be a non-empty string or list of strings.")
        chrom.append(str(record["chrom"]))
        pos.append(_int_option(record["pos"], "pos", minimum=1, maximum=2_000_000_000))
        ref.append(str(record["ref"]))
        ids.append(str(record.get("id") or "."))
        qual.append(float(record["qual"]) if record.get("qual") is not None else None)
        names = record.get("filters") or []
        if isinstance(names, str) or not all(isinstance(item, str) and item.isidentifier() for item in names):
            raise ConfigurationError("filters must be a list of filter names (letters, digits, underscore).")
        failed.append(list(names))
        alts.append(alt_list)
    width = max(len(item) for item in alts)
    alt_array = np.array([item + [""] * (width - len(item)) for item in alts], dtype="S")
    callset = {
        "variants/CHROM": np.array(chrom, dtype="S"),
        "variants/POS": np.array(pos, dtype="i4"),
        "variants/ID": np.array(ids, dtype="S"),
        "variants/REF": np.array(ref, dtype="S"),
        "variants/ALT": alt_array,
        "variants/QUAL": np.array([f"{q:g}" if q is not None else "." for q in qual], dtype="S"),
        "variants/FILTER_PASS": np.array([not item for item in failed], dtype=bool),
    }
    for name in sorted({item for names in failed for item in names}):
        callset[f"variants/FILTER_{name}"] = np.array([name in names for names in failed], dtype=bool)
    allel.write_vcf(str(target), callset, fill={"ALT": b""})
    return _result("write_sites_vcf", output_vcf=str(target), variant_count=len(chrom))


def haploidify_genotypes(
    genotypes: Sequence[Sequence[Sequence[int]]],
    *,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Make pseudo-haplotypes by picking one random allele per sample and variant (``GenotypeArray.haploidify_samples``).

    ``seed`` makes the random choice reproducible by seeding NumPy's legacy generator before the call.
    Returns a ``variants x samples`` list.
    """

    np = _numpy()
    array = _genotypes(genotypes)
    if seed is not None:
        np.random.seed(_int_option(seed, "seed", maximum=2**32 - 1))
    return _result("haploidify_genotypes", seed=seed, haplotypes=_plain(_numpy().asarray(array.haploidify_samples().values)))


def fill_masked_genotypes(
    genotypes: Sequence[Sequence[Sequence[int]]],
    mask: Sequence[Sequence[bool]],
    *,
    fill_value: int = -1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Replace the calls flagged by a ``variants x samples`` boolean mask with ``fill_value`` (``Genotypes.fill_masked``)."""

    np = _numpy()
    allel = _allel()
    array = _int_array(genotypes, "genotypes", ndim=3)
    flags = np.asarray(mask)
    if flags.dtype != bool or flags.shape != array.shape[:2]:
        raise ConfigurationError("mask must be a boolean variants x samples array.")
    fill = _int_option(fill_value, "fill_value", minimum=-1)
    genotype_array = allel.GenotypeArray(array)
    genotype_array.mask = flags
    filled = genotype_array.fill_masked(value=fill)
    return _result("fill_masked_genotypes", fill_value=fill, genotypes=_plain(np.asarray(filled)))


_GENOTYPE_KINDS = ("called", "missing", "het", "hom", "hom_ref", "hom_alt", "call")


def locate_genotype_calls(
    genotypes: Sequence[Sequence[Sequence[int]]],
    *,
    kind: str = "called",
    allele: int | None = None,
    call: Sequence[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Boolean ``variants x samples`` mask of genotype calls of one kind.

    ``kind`` selects the scikit-allel method: ``called`` (is_called), ``missing`` (is_missing), ``het`` (is_het),
    ``hom`` (is_hom; with ``allele`` only homozygous for that allele), ``hom_ref`` (is_hom_ref), ``hom_alt``
    (is_hom_alt) or ``call`` (is_call; ``call`` is the genotype such as ``[0, 1]`` to match).
    """

    array = _genotypes(genotypes)
    chosen = _choice(kind, "kind", _GENOTYPE_KINDS)
    if chosen == "hom":
        allele_value = None if allele is None else _int_option(allele, "allele")
        mask = array.is_hom(allele=allele_value)
    elif chosen == "call":
        if call is None or len(call) != array.shape[2]:
            raise ConfigurationError("call must be a genotype with one allele per ploidy.")
        mask = array.is_call(tuple(_int_option(item, "call allele", minimum=-1) for item in call))
    else:
        mask = getattr(array, f"is_{chosen}")()
    return _result("locate_genotype_calls", kind=chosen, mask=_plain(mask), count=int(mask.sum()))


_HAPLOTYPE_KINDS = ("called", "missing", "ref", "alt", "call")


def locate_haplotype_calls(
    haplotypes: Sequence[Sequence[int]],
    *,
    kind: str = "called",
    allele: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Boolean ``variants x haplotypes`` mask of haplotype calls of one kind.

    ``kind``: ``called`` (is_called), ``missing`` (is_missing), ``ref`` (is_ref), ``alt`` (is_alt) or ``call``
    (is_call; ``allele`` is the allele index to match).
    """

    array = _haplotypes(haplotypes)
    chosen = _choice(kind, "kind", _HAPLOTYPE_KINDS)
    if chosen == "call":
        mask = array.is_call(_int_option(allele, "allele"))
    else:
        mask = getattr(array, f"is_{chosen}")()
    return _result("locate_haplotype_calls", kind=chosen, mask=_plain(mask), count=int(mask.sum()))


_COUNT_GENOTYPE_KINDS = ("hom", "hom_ref", "hom_alt", "het", "call")


def count_genotype_calls(
    genotypes: Sequence[Sequence[Sequence[int]]],
    *,
    kind: str = "hom",
    axis: int | None = None,
    allele: int | None = None,
    call: Sequence[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count genotype calls of one kind in total (``axis`` None), per variant (1) or per sample (0).

    ``kind``: ``hom`` (count_hom; with ``allele``), ``hom_ref``, ``hom_alt``, ``het`` or ``call``
    (count_call; ``call`` is the genotype to match). Note: scikit-allel axis 0 sums over variants, giving one
    value per sample; axis 1 gives one value per variant.
    """

    array = _genotypes(genotypes)
    chosen = _choice(kind, "kind", _COUNT_GENOTYPE_KINDS)
    if axis not in (None, 0, 1) or isinstance(axis, bool):
        raise ConfigurationError("axis must be None, 0 or 1.")
    if chosen == "hom":
        allele_value = None if allele is None else _int_option(allele, "allele")
        counts = array.count_hom(allele=allele_value, axis=axis)
    elif chosen == "call":
        if call is None or len(call) != array.shape[2]:
            raise ConfigurationError("call must be a genotype with one allele per ploidy.")
        counts = array.count_call(tuple(_int_option(item, "call allele", minimum=-1) for item in call), axis=axis)
    else:
        counts = getattr(array, f"count_{chosen}")(axis=axis)
    return _result("count_genotype_calls", kind=chosen, axis=axis, counts=_plain(counts))


def count_haplotype_calls(
    haplotypes: Sequence[Sequence[int]],
    *,
    kind: str = "ref",
    axis: int | None = None,
    allele: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Count haplotype calls of one kind (``ref`` count_ref, ``alt`` count_alt, ``call`` count_call with ``allele``).

    ``axis`` None gives the total, 0 one value per haplotype, 1 one value per variant.
    """

    array = _haplotypes(haplotypes)
    chosen = _choice(kind, "kind", ("ref", "alt", "call"))
    if axis not in (None, 0, 1) or isinstance(axis, bool):
        raise ConfigurationError("axis must be None, 0 or 1.")
    if chosen == "call":
        counts = array.count_call(_int_option(allele, "allele"), axis=axis)
    else:
        counts = getattr(array, f"count_{chosen}")(axis=axis)
    return _result("count_haplotype_calls", kind=chosen, axis=axis, counts=_plain(counts))


def genotype_dosage(
    genotypes: Sequence[Sequence[Sequence[int]]],
    *,
    allele: str = "alt",
    fill: int = 0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Genotypes to the number of reference (``allele='ref'``, to_n_ref) or alternate (``'alt'``, to_n_alt) alleles.

    ``fill`` is the value for missing calls.
    """

    array = _genotypes(genotypes)
    chosen = _choice(allele, "allele", ("ref", "alt"))
    filler = _int_option(fill, "fill", minimum=-1000, maximum=1000)
    dosage = getattr(array, f"to_n_{chosen}")(fill=filler)
    return _result("genotype_dosage", allele=chosen, fill=filler, dosage=_plain(dosage))


def genotype_allele_counts_per_call(
    genotypes: Sequence[Sequence[Sequence[int]]],
    *,
    max_allele: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-call allele counts: for each call, how many times each allele ``0..max_allele`` occurs (``to_allele_counts``)."""

    array = _genotypes(genotypes)
    limit = None if max_allele is None else _int_option(max_allele, "max_allele")
    counts = array.to_allele_counts(max_allele=limit)
    return _result(
        "genotype_allele_counts_per_call",
        max_allele=counts.shape[2] - 1,
        counts=_plain(_numpy().asarray(counts)),
    )


def genotypes_to_gt_strings(
    genotypes: Sequence[Sequence[Sequence[int]]],
    *,
    max_allele: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Format genotype calls as VCF GT strings such as ``0/1`` and ``./.`` (``Genotypes.to_gt``)."""

    array = _genotypes(genotypes)
    limit = None if max_allele is None else _int_option(max_allele, "max_allele")
    gt = array.to_gt(max_allele=limit)
    return _result("genotypes_to_gt_strings", gt=_plain(_numpy().asarray(gt)))


def remap_alleles(
    values: Sequence[Any],
    mapping: Sequence[Sequence[int]],
    *,
    kind: str = "genotypes",
    max_allele: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Recode allele indices through a per-variant ``variants x old_alleles`` mapping (``map_alleles``).

    ``kind``: ``genotypes`` (variants x samples x ploidy), ``haplotypes`` (variants x haplotypes) or
    ``allele_counts`` (variants x alleles; ``max_allele`` limits the new width). A mapped value of -1 drops the allele.
    """

    chosen = _choice(kind, "kind", ("genotypes", "haplotypes", "allele_counts"))
    table = _int_array(mapping, "mapping", ndim=2)
    if chosen == "genotypes":
        array = _genotypes(values, "values")
    elif chosen == "haplotypes":
        array = _haplotypes(values, "values")
    else:
        array = _allele_counts(values, "values")
    if table.shape[0] != array.shape[0]:
        raise ConfigurationError("mapping needs one row per variant.")
    if chosen == "allele_counts":
        limit = None if max_allele is None else _int_option(max_allele, "max_allele")
        out = array.map_alleles(table, max_allele=limit)
    else:
        out = array.map_alleles(table)
    return _result("remap_alleles", kind=chosen, values=_plain(_numpy().asarray(out)))


def summarize_alleles_per_variant(
    allele_counts: Sequence[Sequence[int]],
    *,
    kind: str = "allelism",
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-variant integer summary of allele counts.

    ``kind='allelism'`` is the number of distinct alleles observed (``allelism``); ``kind='max_allele'`` is the
    highest allele index observed (``max_allele``, -1 when nothing is observed).
    """

    array = _allele_counts(allele_counts)
    chosen = _choice(kind, "kind", ("allelism", "max_allele"))
    return _result("summarize_alleles_per_variant", kind=chosen, values=_plain(getattr(array, chosen)()))


_SITE_KINDS = (
    "variant",
    "non_variant",
    "segregating",
    "non_segregating",
    "singleton",
    "doubleton",
    "biallelic",
    "biallelic_01",
)


def locate_allele_count_sites(
    allele_counts: Sequence[Sequence[int]],
    *,
    kind: str = "variant",
    allele: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Boolean per-variant mask of sites of one kind from allele counts.

    ``kind``: ``variant`` (is_variant), ``non_variant``, ``segregating``, ``non_segregating``, ``singleton`` and
    ``doubleton`` (with ``allele``), ``biallelic`` or ``biallelic_01`` (biallelic for alleles 0 and 1).
    """

    array = _allele_counts(allele_counts)
    chosen = _choice(kind, "kind", _SITE_KINDS)
    if chosen in ("singleton", "doubleton"):
        mask = getattr(array, f"is_{chosen}")(_int_option(allele, "allele"))
    else:
        mask = getattr(array, f"is_{chosen}")()
    return _result("locate_allele_count_sites", kind=chosen, mask=_plain(mask), count=int(mask.sum()))


def count_allele_count_sites(
    allele_counts: Sequence[Sequence[int]],
    *,
    kind: str = "variant",
    allele: int = 1,
    verbose: bool = True,
) -> dict[str, Any]:
    """Number of sites of one kind: ``variant``, ``non_variant``, ``segregating``, ``non_segregating``, ``singleton``, ``doubleton``.

    ``allele`` applies to ``singleton`` and ``doubleton`` (count_singleton, count_doubleton).
    """

    array = _allele_counts(allele_counts)
    chosen = _choice(kind, "kind", _SITE_KINDS[:6])
    if chosen in ("singleton", "doubleton"):
        total = getattr(array, f"count_{chosen}")(_int_option(allele, "allele"))
    else:
        total = getattr(array, f"count_{chosen}")()
    return _result("count_allele_count_sites", kind=chosen, count=int(total))


def feature_coverage_mask(
    features: Sequence[Sequence[int]],
    size: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Boolean mask along a contig of the positions covered by features (``FeatureTable.to_mask``).

    ``features`` is a list of 1-based inclusive ``[start, end]`` pairs and ``size`` the contig length; the mask has
    ``size + 1`` entries (index 0 unused, as in scikit-allel).
    """

    np = _numpy()
    allel = _allel()
    length = _int_option(size, "size", minimum=1, maximum=500_000_000)
    pairs = _int_array(features, "features", ndim=2, minimum=0, maximum=_MAX_POSITION)
    if pairs.shape[1] != 2 or (pairs[:, 1] < pairs[:, 0]).any():
        raise ConfigurationError("features must be [start, end] pairs with end >= start.")
    table = allel.FeatureTable(
        np.rec.fromarrays([pairs[:, 0], pairs[:, 1]], names="start,end"),
    )
    mask = table.to_mask(length)
    return _result("feature_coverage_mask", size=length, mask=_plain(mask), covered=int(mask.sum()))


def locate_sorted_positions(
    positions: Sequence[int],
    keys: int | Sequence[int],
    *,
    strict: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Locate one position (``locate_key``) or several (``locate_keys``) in a sorted position index.

    A single integer key returns its index; a list returns a boolean mask over ``positions``. With ``strict=False`` missing keys
    are ignored instead of raising. Positions must be unique and sorted.
    """

    allel = _allel()
    index = allel.SortedIndex(_positions(positions, "positions"))
    if not index.is_unique:
        raise ConfigurationError("positions must be unique.")
    try:
        if isinstance(keys, int) and not isinstance(keys, bool):
            return _result("locate_sorted_positions", mode="key", located=int(index.locate_key(keys)))
        key_array = _int_array(keys, "keys", ndim=1, minimum=1, maximum=_MAX_POSITION)
        mask = index.locate_keys(key_array, strict=bool(strict))
    except KeyError as exc:
        raise ConfigurationError(f"Position not found in the index: {exc}") from exc
    return _result("locate_sorted_positions", mode="keys", mask=_plain(mask), count=int(mask.sum()))


def intersect_sorted_positions(
    positions: Sequence[int],
    other: Sequence[int],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Intersect two sets of sorted genomic positions (``locate_intersection`` masks and ``intersect`` values)."""

    allel = _allel()
    first = allel.SortedIndex(_positions(positions, "positions"))
    second = allel.SortedIndex(_positions(other, "other"))
    if not (first.is_unique and second.is_unique):
        raise ConfigurationError("positions and other must each be unique.")
    first_mask, second_mask = first.locate_intersection(second)
    return _result(
        "intersect_sorted_positions",
        positions=_plain(first.intersect(second).values),
        mask_positions=_plain(first_mask),
        mask_other=_plain(second_mask),
    )


def _range_bounds(ranges: object) -> tuple[Any, Any]:
    pairs = _int_array(ranges, "ranges", ndim=2, minimum=1, maximum=_MAX_POSITION)
    if pairs.shape[1] != 2 or (pairs[:, 1] < pairs[:, 0]).any():
        raise ConfigurationError("ranges must be [start, stop] pairs (1-based, inclusive) with stop >= start.")
    return pairs[:, 0], pairs[:, 1]


def select_positions_in_ranges(
    positions: Sequence[int],
    ranges: Sequence[Sequence[int]],
    *,
    strict: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Positions that fall inside one or several 1-based inclusive ranges (``locate_range``, ``intersect_range``, ``locate_ranges``, ``intersect_ranges``).

    ``ranges`` is a list of ``[start, stop]`` pairs; one pair equals the single-range methods. Returns the boolean mask over
    ``positions`` and the selected positions. With ``strict=True`` a range without any position raises.
    """

    allel = _allel()
    index = allel.SortedIndex(_positions(positions, "positions"))
    starts, stops = _range_bounds(ranges)
    try:
        mask = index.locate_ranges(starts, stops, strict=bool(strict))
    except KeyError as exc:
        raise ConfigurationError(f"A range contains no position: {exc}") from exc
    return _result(
        "select_positions_in_ranges",
        mask=_plain(mask),
        positions=_plain(index.values[mask]),
        count=int(mask.sum()),
    )


def overlap_positions_with_ranges(
    positions: Sequence[int],
    ranges: Sequence[Sequence[int]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mask of positions inside any range and mask of ranges that contain at least one position (``locate_intersection_ranges``)."""

    allel = _allel()
    index = allel.SortedIndex(_positions(positions, "positions"))
    starts, stops = _range_bounds(ranges)
    position_mask, range_mask = index.locate_intersection_ranges(starts, stops)
    return _result(
        "overlap_positions_with_ranges",
        position_mask=_plain(position_mask),
        range_mask=_plain(range_mask),
    )


def create_allele_index_mapping(
    ref: Sequence[str],
    alt: Sequence[Sequence[str]],
    alleles: Sequence[Sequence[str]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mapping from VCF allele indices (REF, ALT...) to indices in a per-variant allele list such as ancestral/derived (``create_allele_mapping``).

    ``ref`` has one base string per variant, ``alt`` one list of alternate alleles per variant (padded
    with ``''``), ``alleles`` the target alleles per variant. Alleles absent from the target list map to -1.
    """

    np = _numpy()
    allel = _allel()
    try:
        ref_array = np.array(ref, dtype="U")
        alt_array = np.array(alt, dtype="U")
        allele_array = np.array(alleles, dtype="U")
    except (TypeError, ValueError) as exc:
        raise ConfigurationError("ref, alt and alleles must be rectangular string arrays.") from exc
    if ref_array.ndim != 1 or alt_array.ndim != 2 or allele_array.ndim != 2:
        raise ConfigurationError("ref must be 1-D, alt and alleles 2-D.")
    if not ref_array.shape[0] == alt_array.shape[0] == allele_array.shape[0]:
        raise ConfigurationError("ref, alt and alleles need one row per variant.")
    mapping = allel.create_allele_mapping(ref_array, alt_array, allele_array)
    return _result("create_allele_index_mapping", mapping=_plain(mapping))


def principal_coordinates_analysis(
    distance_matrix: Sequence[Sequence[float]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Principal coordinate analysis of a square genetic distance matrix (``allel.pcoa``); returns coordinates and explained variance ratios."""

    np = _numpy()
    allel = _allel()
    try:
        matrix = np.asarray(distance_matrix, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError("distance_matrix must be numeric.") from exc
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] < 3 or matrix.size > 4_000_000:
        raise ConfigurationError("distance_matrix must be a square matrix with at least 3 rows.")
    if not np.isfinite(matrix).all() or not np.allclose(matrix, matrix.T):
        raise ConfigurationError("distance_matrix must be finite and symmetric.")
    coords, explained = allel.pcoa(matrix)
    return _result("principal_coordinates_analysis", coordinates=_plain(coords), explained_ratio=_plain(explained))


def _states(x: object, states: object, pos: object) -> tuple[Any, Any, Any]:
    np = _numpy()
    values = _int_array(x, "x", ndim=1, minimum=-1000, maximum=_MAX_COUNT)
    chosen = _int_array(states, "states", ndim=1, minimum=-1000, maximum=_MAX_COUNT)
    positions = None
    if pos is not None:
        positions = _int_array(pos, "pos", ndim=1, minimum=1, maximum=_MAX_POSITION)
        if positions.shape != values.shape or (np.diff(positions) < 0).any():
            raise ConfigurationError("pos must be sorted and as long as x.")
    return values, chosen, positions


def _table(frame: Any) -> list[dict[str, Any]]:
    return [{str(k): _plain(v) for k, v in row.items()} for row in frame.to_dict(orient="records")]


def tabulate_state_transitions(
    x: Sequence[int],
    states: Sequence[int],
    *,
    pos: Sequence[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Table of transitions between the listed states along a series (``allel.tabulate_state_transitions``).

    ``x`` is the state series, ``states`` the states of interest, ``pos`` optional sorted 1-based positions.
    Rows hold lstate, rstate, lidx, ridx and, with ``pos``, lpos and rpos.
    """

    values, chosen, positions = _states(x, states, pos)
    frame = _allel().tabulate_state_transitions(values, chosen, pos=positions)
    return _result("tabulate_state_transitions", transitions=_table(frame))


def tabulate_state_blocks(
    x: Sequence[int],
    states: Sequence[int],
    *,
    pos: Sequence[int] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Table of contiguous blocks of the listed states along a series (``allel.tabulate_state_blocks``).

    Rows hold the state, support, start/stop indices, positions with ``pos``, and the minimum and maximum
    block length. -1 marks an undefined boundary.
    """

    values, chosen, positions = _states(x, states, pos)
    frame = _allel().tabulate_state_blocks(values, chosen, pos=positions)
    return _result("tabulate_state_blocks", blocks=_table(frame))


def scale_genotype_dosages(
    dosages: Sequence[Sequence[float]],
    *,
    method: str = "patterson",
    ploidy: int = 2,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scale a ``variants x samples`` alternate-allele dosage matrix per variant before PCA.

    ``method='center'`` subtracts the variant mean (``CenterScaler``); ``'patterson'`` also divides by
    sqrt(p(1-p)) with p the allele frequency (``PattersonScaler`` with ``ploidy``). Monomorphic variants give undefined values, returned as null.
    """

    np = _numpy()
    allel = _allel()
    chosen = _choice(method, "method", ("center", "patterson"))
    ploid = _int_option(ploidy, "ploidy", minimum=1, maximum=16)
    try:
        matrix = np.asarray(dosages, dtype=float)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError("dosages must be numeric.") from exc
    if matrix.ndim != 2 or matrix.size == 0 or matrix.size > _MAX_ITEMS or not np.isfinite(matrix).all():
        raise ConfigurationError("dosages must be a finite, non-empty 2-D array.")
    scaler = allel.CenterScaler() if chosen == "center" else allel.PattersonScaler(ploidy=ploid)
    scaled = scaler.fit(matrix).transform(matrix)
    return _result("scale_genotype_dosages", method=chosen, ploidy=ploid, scaled=_plain(scaled))


def fold_site_frequency_spectrum(
    spectrum: Sequence[int],
    sample_size: int,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Fold a one-population site frequency spectrum onto minor-allele counts (``allel.fold_sfs``); ``sample_size`` is the chromosome count n and a shorter spectrum is zero-padded to n + 1 entries. As in scikit-allel, bin k is s[k] + s[n - k] for k < (n + 1) // 2, so use an odd n for a symmetric fold."""

    np = _numpy()
    allel = _allel()
    counts = _int_array(spectrum, "spectrum", ndim=1, minimum=0, maximum=_MAX_COUNT)
    n = _int_option(sample_size, "sample_size", minimum=1, maximum=100_000)
    if counts.shape[0] > n + 1:
        raise ConfigurationError("spectrum must have at most sample_size + 1 entries.")
    folded = allel.fold_sfs(counts, n)
    return _result("fold_site_frequency_spectrum", sample_size=n, folded=_plain(np.asarray(folded)))


__all__ = [
    "count_allele_count_sites",
    "count_genotype_calls",
    "count_haplotype_calls",
    "create_allele_index_mapping",
    "feature_coverage_mask",
    "fill_masked_genotypes",
    "fold_site_frequency_spectrum",
    "genotype_allele_counts_per_call",
    "genotype_dosage",
    "genotypes_to_gt_strings",
    "haploidify_genotypes",
    "intersect_sorted_positions",
    "locate_allele_count_sites",
    "locate_genotype_calls",
    "locate_haplotype_calls",
    "locate_sorted_positions",
    "overlap_positions_with_ranges",
    "principal_coordinates_analysis",
    "remap_alleles",
    "scale_genotype_dosages",
    "select_positions_in_ranges",
    "summarize_alleles_per_variant",
    "tabulate_state_blocks",
    "tabulate_state_transitions",
    "write_sites_vcf",
]

_register_public_verbose(globals(), tuple(__all__))
