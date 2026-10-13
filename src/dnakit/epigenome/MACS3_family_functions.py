"""MACS3 wrappers for pileup and variant-model variants not reached by the existing wrappers.

The functions call MACS3's Cython functions directly. They cover count-weighted paired-end
pileups, pileups centred on fragment ends, the HMMRATAC length-weighted pileup, and the
allele-specific heterozygous model with its genotype quality. MACS3 is imported only when a
function runs. A missing package raises ``BackendUnavailableError`` with the
``pip install "dnakit[bio-macs3]"`` hint.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.epigenome.MACS3_functions import (
    _MACS3_EXTRA,
    _count_array,
    _integer,
    _number,
    _pileup_result,
    _plain,
    _positions,
    _provenance,
    _require_sorted,
)
from dnakit.exceptions import ConfigurationError

_MAX_COUNT = 65_535
_MAX_WINDOW = 1_000_000


def _fragments(left: object, right: object) -> tuple[list[int], list[int]]:
    starts = _positions(left, "left_positions")
    ends = _positions(right, "right_positions")
    if len(starts) != len(ends) or not starts:
        raise ConfigurationError(
            "left_positions and right_positions must have equal, non-zero length."
        )
    if any(end <= start for start, end in zip(starts, ends, strict=True)):
        raise ConfigurationError(
            "Every fragment must have its right coordinate after its left one."
        )
    return starts, ends


def _counts(value: object, expected: int) -> list[int]:
    items = _positions(value, "counts")
    if len(items) != expected:
        raise ConfigurationError("counts must have one value per fragment.")
    if any(item > _MAX_COUNT for item in items):
        raise ConfigurationError(f"counts must be at most {_MAX_COUNT}.")
    return items


def _lrc_array(numpy: Any, starts: list[int], ends: list[int], counts: list[int]) -> Any:
    return numpy.array(
        list(zip(starts, ends, counts, strict=True)),
        dtype=[("l", "<i4"), ("r", "<i4"), ("c", "<u2")],
    )


def build_counted_fragment_pileup(
    left_positions: Sequence[int],
    right_positions: Sequence[int],
    counts: Sequence[int],
    *,
    scale_factor: float = 1.0,
    baseline_value: float = 0.0,
    left_sorted: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pile up paired-end fragments that each stand for several identical reads.

    Calls ``MACS3.Signal.PileupV2.pileup_from_LRC_as_list``: fragment ``i`` covers
    ``[left_positions[i], right_positions[i])`` and adds ``counts[i]`` to every base it covers.
    ``left_sorted=True`` promises that the fragments are already sorted by left position (it is
    checked) and lets MACS3 skip its sort. The result has the same layout as
    ``build_paired_fragment_pileup``: ``positions`` are segment ends and ``values`` the pileup of
    the segment that ends at each position.
    """

    starts, ends = _fragments(left_positions, right_positions)
    weights = _counts(counts, len(starts))
    scale = _number(scale_factor, "scale_factor", 0.0, 1_000_000.0)
    baseline = _number(baseline_value, "baseline_value", -1_000_000.0, 1_000_000.0)
    if not isinstance(left_sorted, bool):
        raise ConfigurationError("left_sorted must be a bool.")
    if left_sorted:
        _require_sorted(starts, "left_positions")
    numpy = require("numpy", _MACS3_EXTRA)
    pileup_module = require("MACS3.Signal.PileupV2", _MACS3_EXTRA)
    array = _lrc_array(numpy, starts, ends, weights)
    pileup = pileup_module.pileup_from_LRC_as_list(array, scale, baseline, left_sorted)
    return {
        "status": "completed",
        "operation": "build_counted_fragment_pileup",
        "fragments": len(starts),
        "total_count": sum(weights),
        **_pileup_result(pileup),
        "provenance": _provenance(
            "MACS3",
            "count-weighted pileup of fragments; positions are segment ends, values cover the "
            "segment that ends at each position",
            source_function="MACS3.Signal.PileupV2.pileup_from_LRC_as_list",
        ),
    }


def build_fragment_end_pileup(
    left_positions: Sequence[int],
    right_positions: Sequence[int],
    counts: Sequence[int],
    *,
    window: int,
    scale_factor: float = 1.0,
    baseline_value: float = 0.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pile up windows centred on both ends of count-weighted fragments.

    Calls ``MACS3.Signal.PileupV2.pileup_from_LRC_centers_as_list``. Every fragment contributes two
    windows of ``window`` bases with weight ``counts[i]``: one starting at
    ``left - window // 2`` and one at ``right - window // 2``. The fragment body is not covered, so
    the pileup shows where fragments start and end. The layout of the result is the same as
    ``build_counted_fragment_pileup``.
    """

    starts, ends = _fragments(left_positions, right_positions)
    weights = _counts(counts, len(starts))
    width = _integer(window, "window", 1, _MAX_WINDOW)
    scale = _number(scale_factor, "scale_factor", 0.0, 1_000_000.0)
    baseline = _number(baseline_value, "baseline_value", -1_000_000.0, 1_000_000.0)
    numpy = require("numpy", _MACS3_EXTRA)
    pileup_module = require("MACS3.Signal.PileupV2", _MACS3_EXTRA)
    array = _lrc_array(numpy, starts, ends, weights)
    pileup = pileup_module.pileup_from_LRC_centers_as_list(array, width, scale, baseline)
    return {
        "status": "completed",
        "operation": "build_fragment_end_pileup",
        "fragments": len(starts),
        "window": width,
        "total_count": sum(weights),
        **_pileup_result(pileup),
        "provenance": _provenance(
            "MACS3",
            "count-weighted pileup of windows centred on the left and right end of each fragment",
            source_function="MACS3.Signal.PileupV2.pileup_from_LRC_centers_as_list",
        ),
    }


def build_hmmratac_fragment_pileup(
    left_positions: Sequence[int],
    right_positions: Sequence[int],
    length_weights: Mapping[str, float],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pile up paired-end fragments weighted by their length, as HMMRATAC does.

    Calls ``MACS3.Signal.PileupV2.pileup_from_LR_hmmratac``. ``length_weights`` maps a fragment
    length (``right - left``, written as a decimal string such as ``"40"``) to the weight added to
    every base of that fragment. Every fragment length present must have a weight. The weights can
    come from ``hmmratac_fragment_weights``. The layout of the result is the same as
    ``build_paired_fragment_pileup``.
    """

    starts, ends = _fragments(left_positions, right_positions)
    if not isinstance(length_weights, Mapping) or not length_weights:
        raise ConfigurationError("length_weights must be a non-empty mapping of length to weight.")
    table: dict[int, float] = {}
    for key, weight in length_weights.items():
        try:
            length = int(key)
        except (TypeError, ValueError) as exc:
            raise ConfigurationError("length_weights keys must be integer lengths.") from exc
        if length < 1:
            raise ConfigurationError("length_weights keys must be positive lengths.")
        table[length] = _number(weight, "length_weights values", 0.0, 1_000_000.0)
    missing = sorted({end - start for start, end in zip(starts, ends, strict=True)} - set(table))
    if missing:
        raise ConfigurationError(f"length_weights has no weight for fragment lengths {missing}.")
    numpy = require("numpy", _MACS3_EXTRA)
    pileup_module = require("MACS3.Signal.PileupV2", _MACS3_EXTRA)
    array = numpy.array(list(zip(starts, ends, strict=True)), dtype=[("l", "<i4"), ("r", "<i4")])
    pileup = pileup_module.pileup_from_LR_hmmratac(array, table)
    positions = [int(_plain(item["p"])) for item in pileup]
    values = [float(_plain(item["v"])) for item in pileup]
    return {
        "status": "completed",
        "operation": "build_hmmratac_fragment_pileup",
        "fragments": len(starts),
        "positions": positions,
        "values": values,
        "provenance": _provenance(
            "MACS3",
            "pileup of fragments weighted by a per-length weight; positions are segment ends, "
            "values cover the segment that ends at each position",
            source_function="MACS3.Signal.PileupV2.pileup_from_LR_hmmratac",
        ),
    }


def calculate_allele_specific_model_likelihood(
    top1_bq_t: Sequence[int],
    top1_bq_c: Sequence[int],
    top2_bq_t: Sequence[int],
    top2_bq_c: Sequence[int],
    *,
    max_allowed_ratio: float = 0.99,
    verbose: bool = True,
) -> dict[str, Any]:
    """Log-likelihood and BIC of the heterozygous model with an allele-specific ratio.

    Calls ``MACS3.Signal.VariantStat.CalModel_Heter_AS``. The treatment reads (``_t``) are fitted
    with a free allele ratio capped at ``max_allowed_ratio``; the control reads (``_c``) are fitted
    without allele preference. Each argument lists the base qualities of the reads that support
    the top allele (``top1``) or the second allele (``top2``). The treatment arrays must hold at
    least one read. ``calculate_allele_model_likelihood`` covers the models without a free ratio.
    """

    ratio = _number(max_allowed_ratio, "max_allowed_ratio", 0.5, 1.0)
    numpy = require("numpy", _MACS3_EXTRA)
    arrays = [
        _count_array(top1_bq_t, "top1_bq_t", numpy),
        _count_array(top1_bq_c, "top1_bq_c", numpy),
        _count_array(top2_bq_t, "top2_bq_t", numpy),
        _count_array(top2_bq_c, "top2_bq_c", numpy),
    ]
    if any(len(array) != len(arrays[0]) for array in arrays):
        raise ConfigurationError("All four base-quality arrays must have the same length.")
    variant = require("MACS3.Signal.VariantStat", _MACS3_EXTRA)
    log_likelihood, bic = variant.CalModel_Heter_AS(*arrays, ratio)
    return {
        "status": "completed",
        "operation": "calculate_allele_specific_model_likelihood",
        "reads": len(arrays[0]),
        "max_allowed_ratio": ratio,
        "log_likelihood": _plain(log_likelihood),
        "bic": _plain(bic),
        "provenance": _provenance(
            "MACS3",
            "heterozygous genotype model with a free allele ratio for the treatment reads",
            source_function="MACS3.Signal.VariantStat.CalModel_Heter_AS",
        ),
    }


def calculate_allele_specific_genotype_quality(
    log_likelihood_1: float,
    log_likelihood_2: float,
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Phred-like score for the allele-specific model against a second model.

    Calls ``MACS3.Signal.VariantStat.calculate_GQ_heterASsig``. With ``L1 = 1`` and
    ``L2 = min(1, exp(log_likelihood_2 - log_likelihood_1))`` (at least ``1e-110``) the score is
    ``int(-10 * log10(L2 / (L1 + L2)))``, or 255 when that ratio is not above ``1e-110``. It is
    3 when the two log-likelihoods are equal and grows as model 1 fits better. The three-model
    ``calculate_genotype_quality`` uses a different formula.
    """

    first = _number(log_likelihood_1, "log_likelihood_1", -1e9, 1e9)
    second = _number(log_likelihood_2, "log_likelihood_2", -1e9, 1e9)
    variant = require("MACS3.Signal.VariantStat", _MACS3_EXTRA)
    score = variant.calculate_GQ_heterASsig(first, second)
    return {
        "status": "completed",
        "operation": "calculate_allele_specific_genotype_quality",
        "log_likelihood_1": first,
        "log_likelihood_2": second,
        "score": int(score),
        "provenance": _provenance(
            "MACS3",
            "int(-10 * log10(L2 / (L1 + L2))) with L1 = 1 and L2 = exp(lnL2 - lnL1) clipped to "
            "[1e-110, 1]",
            source_function="MACS3.Signal.VariantStat.calculate_GQ_heterASsig",
        ),
    }


__all__ = [
    "build_counted_fragment_pileup",
    "build_fragment_end_pileup",
    "build_hmmratac_fragment_pileup",
    "calculate_allele_specific_genotype_quality",
    "calculate_allele_specific_model_likelihood",
]

_register_public_verbose(globals(), tuple(__all__))
