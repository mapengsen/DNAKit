"""Biopython ``Bio.motifs`` operations added by integration pack P19.

A motif is given either as a list of equal-length DNA ``instances`` (binding sites) or as ``counts``, a
mapping of base (A, C, G, T) to a list of per-position counts. ``pseudocounts`` (a number added to every
cell, or a mapping per base) and ``background`` (base probabilities, uniform by default) are applied as in
``Bio.motifs.Motif``. Biopython's defaults are kept: pseudocounts are zero, so a base never seen at a
position scores minus infinity (reported as ``null`` in JSON results). Biopython is imported only when a
function runs.

``weblogo_service_logo`` calls the public WebLogo web service (https://weblogo.threeplusone.com) and needs
network access.
"""

from __future__ import annotations

import importlib.metadata
import math
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-formats"
_PACKAGE = "biopython"
_BASES = "ACGT"
_MAX_WIDTH = 200
_MAX_INSTANCES = 100_000
_MAX_SEQUENCE = 50_000_000
_MAX_DISTRIBUTION_POINTS = 400_000
_THRESHOLD_METHODS = ("fpr", "fnr", "balanced", "patser")
_LOGO_FORMATS = ("png", "png_print", "pdf", "jpeg", "eps", "svg", "txt")
_LOGO_OPTIONS = {
    "stack_width": str,
    "stacks_per_line": str,
    "unit_name": str,
    "first_index": str,
    "logo_start": str,
    "logo_end": str,
    "composition": str,
    "percentCG": str,
    "logo_title": str,
    "logo_label": str,
    "xaxis_label": str,
    "yaxis_label": str,
    "yaxis_scale": str,
    "yaxis_tic_interval": str,
    "color_scheme": str,
    "scale_width": bool,
    "show_errorbars": bool,
    "show_xaxis": bool,
    "show_yaxis": bool,
    "show_ends": bool,
    "show_fineprint": bool,
    "ignore_lower_case": bool,
}


def _provenance(**notes: str) -> dict[str, Any]:
    try:
        installed: str | None = importlib.metadata.version(_PACKAGE)
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": _PACKAGE, "package_version": installed, **notes}


def _bio(module: str) -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return require(module, _EXTRA)


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _finite(value: float) -> float | None:
    """JSON-safe score: infinite or NaN values become None."""

    number = float(value)
    return number if math.isfinite(number) else None


def _row(values: Sequence[float]) -> list[float | None]:
    return [_finite(item) for item in values]


def _table(matrix: Any) -> dict[str, list[float | None]]:
    return {base: _row(list(matrix[base])) for base in _BASES}


def _base_map(value: object, name: str, *, positive_sum: bool = False) -> dict[str, float] | None:
    if value is None:
        return None
    if isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a number or a mapping of base to number.")
    if isinstance(value, (int, float)):
        if not math.isfinite(value) or value < 0:
            raise ConfigurationError(f"{name} must be a non-negative number.")
        return dict.fromkeys(_BASES, float(value))
    if not isinstance(value, Mapping) or set(value) != set(_BASES):
        raise ConfigurationError(f"{name} must give a value for each of A, C, G and T.")
    result = {}
    for base in _BASES:
        item = value[base]
        if isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(item) or item < 0:
            raise ConfigurationError(f"{name}[{base!r}] must be a non-negative number.")
        result[base] = float(item)
    if positive_sum and sum(result.values()) <= 0:
        raise ConfigurationError(f"{name} must have a positive sum.")
    return result


def _motif(
    instances: Sequence[str] | None,
    counts: Mapping[str, Sequence[float]] | None,
    pseudocounts: Any,
    background: Any,
) -> Any:
    if (instances is None) == (counts is None):
        raise ConfigurationError("Provide exactly one of instances and counts.")
    motifs = _bio("Bio.motifs")
    seq = _bio("Bio.Seq")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if instances is not None:
                if isinstance(instances, (str, bytes)) or not isinstance(instances, Sequence) or not instances:
                    raise ConfigurationError("instances must be a non-empty list of DNA strings.")
                if len(instances) > _MAX_INSTANCES:
                    raise ConfigurationError(f"instances may hold at most {_MAX_INSTANCES} sites.")
                cleaned = []
                for site in instances:
                    if not isinstance(site, str) or not site or set(site.upper()) - set(_BASES):
                        raise ConfigurationError("instances must be DNA strings of A, C, G and T.")
                    cleaned.append(seq.Seq(site.upper()))
                if len({len(site) for site in cleaned}) != 1 or len(cleaned[0]) > _MAX_WIDTH:
                    raise ConfigurationError(
                        f"instances must all have the same length (at most {_MAX_WIDTH})."
                    )
                motif = motifs.create(cleaned)
            else:
                if not isinstance(counts, Mapping) or set(counts) != set(_BASES):
                    raise ConfigurationError("counts must give a list for each of A, C, G and T.")
                widths = {len(counts[base]) for base in _BASES}
                if len(widths) != 1 or not 1 <= next(iter(widths)) <= _MAX_WIDTH:
                    raise ConfigurationError(f"counts lists must have the same length (1 to {_MAX_WIDTH}).")
                table = {}
                for base in _BASES:
                    table[base] = []
                    for item in counts[base]:
                        if isinstance(item, bool) or not isinstance(item, (int, float)) or not math.isfinite(item) or item < 0:
                            raise ConfigurationError("counts must be non-negative numbers.")
                        table[base].append(float(item))
                motif = motifs.Motif(alphabet=_BASES, counts=table)
            chosen_background = _base_map(background, "background", positive_sum=True)
            if chosen_background is not None:
                total = sum(chosen_background.values())
                motif.background = {base: value / total for base, value in chosen_background.items()}
            chosen_pseudocounts = _base_map(pseudocounts, "pseudocounts")
            if chosen_pseudocounts is not None:
                motif.pseudocounts = chosen_pseudocounts
    except (ValueError, KeyError) as exc:
        raise ConfigurationError(f"Biopython could not build the motif: {exc}") from exc
    return motif


def _pssm(motif: Any) -> Any:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return motif.pssm


def motif_matrices(
    instances: Sequence[str] | None = None,
    *,
    counts: Mapping[str, Sequence[float]] | None = None,
    pseudocounts: float | Mapping[str, float] | None = None,
    background: Mapping[str, float] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a motif from binding sites and return its count, probability and log-odds matrices.

    ``Bio.motifs.create`` turns the ``instances`` (or the given ``counts``) into a count matrix;
    ``Motif.pwm`` normalises it into a position weight matrix (probabilities, with ``pseudocounts``) and
    ``Motif.pssm`` gives the log2-odds position-specific scoring matrix against ``background`` (uniform by
    default). Also returns the consensus, the degenerate (IUPAC) consensus and the anticonsensus.
    Matrices are ``{base: [value per position]}``; ``null`` marks minus infinity.
    """

    motif = _motif(instances, counts, pseudocounts, background)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pwm = motif.pwm
        pssm = motif.pssm
        result = {
            "status": "completed",
            "operation": "motif_matrices",
            "width": int(motif.length),
            "instance_count": len(instances) if instances is not None else None,
            "counts": _table(motif.counts),
            "pwm": _table(pwm),
            "pssm": _table(pssm),
            "background": {base: float(motif.background[base]) for base in _BASES},
            "consensus": str(motif.consensus),
            "degenerate_consensus": str(motif.degenerate_consensus),
            "anticonsensus": str(motif.anticonsensus),
        }
    result["provenance"] = _provenance(
        parser="Bio.motifs.create; Bio.motifs.Motif.pwm; Bio.motifs.Motif.pssm",
        infinity="null marks -inf (a base with probability 0 and no pseudocount)",
    )
    return result


def motif_anticonsensus(
    instances: Sequence[str] | None = None,
    *,
    counts: Mapping[str, Sequence[float]] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the anticonsensus of a motif: the least frequent base at every position.

    This is ``Bio.motifs.Motif.anticonsensus``. It is the sequence that minimises the match to the motif;
    it is also returned together with the consensus and the degenerate consensus for comparison.
    """

    motif = _motif(instances, counts, None, None)
    return {
        "status": "completed",
        "operation": "motif_anticonsensus",
        "anticonsensus": str(motif.anticonsensus),
        "consensus": str(motif.consensus),
        "degenerate_consensus": str(motif.degenerate_consensus),
        "provenance": _provenance(parser="Bio.motifs.Motif.anticonsensus"),
    }


def pssm_search(
    sequence: str,
    instances: Sequence[str] | None = None,
    *,
    counts: Mapping[str, Sequence[float]] | None = None,
    threshold: float = 0.0,
    both_strands: bool = True,
    pseudocounts: float | Mapping[str, float] | None = 0.5,
    background: Mapping[str, float] | None = None,
    max_hits: int = 100_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan a DNA sequence with a motif's log-odds PSSM and report hits at or above ``threshold``.

    This is ``Bio.motifs.matrix.PositionSpecificScoringMatrix.search``. The sequence is upper-cased as
    Biopython does. Each hit has ``position`` (Biopython's convention: counted from the start on the
    forward strand, negative for the reverse strand, where ``position + len(sequence)`` is the forward
    start), the forward-strand 0-based half-open ``start``/``end``, the ``strand`` and the ``score``.
    ``pseudocounts`` defaults to 0.5 here so that unseen bases do not give -inf; pass 0 for the
    Biopython default.
    """

    if not isinstance(sequence, str) or not sequence or len(sequence) > _MAX_SEQUENCE:
        raise ConfigurationError("sequence must be a non-empty DNA string.")
    if set(sequence.upper()) - set("ACGTN"):
        raise ConfigurationError("sequence may contain only A, C, G, T and N.")
    if isinstance(threshold, bool) or not isinstance(threshold, (int, float)) or not math.isfinite(threshold):
        raise ConfigurationError("threshold must be a finite number.")
    limit = _integer(max_hits, "max_hits", minimum=1, maximum=10_000_000)
    motif = _motif(instances, counts, pseudocounts, background)
    width = int(motif.length)
    if width > len(sequence):
        raise ConfigurationError("the motif is longer than the sequence.")
    seq = _bio("Bio.Seq")
    pssm = _pssm(motif)
    hits = []
    truncated = False
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for position, score in pssm.search(
            seq.Seq(sequence), threshold=float(threshold), both=_flag(both_strands, "both_strands")
        ):
            if len(hits) >= limit:
                truncated = True
                break
            position = int(position)
            reverse = position < 0
            start = position + len(sequence) if reverse else position
            hits.append(
                {
                    "position": position,
                    "start": start,
                    "end": start + width,
                    "strand": "-" if reverse else "+",
                    "score": float(score),
                }
            )
    return {
        "status": "completed",
        "operation": "pssm_search",
        "width": width,
        "threshold": float(threshold),
        "hit_count": len(hits),
        "truncated": truncated,
        "hits": hits,
        "provenance": _provenance(parser="Bio.motifs.matrix.PositionSpecificScoringMatrix.search"),
    }


def pssm_score_statistics(
    instances: Sequence[str] | None = None,
    *,
    counts: Mapping[str, Sequence[float]] | None = None,
    pseudocounts: float | Mapping[str, float] | None = 0.5,
    background: Mapping[str, float] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Mean, standard deviation, maximum and minimum of a motif's PSSM score under the background.

    ``mean`` and ``std`` are the expected value and standard deviation of the score of a random sequence
    drawn from ``background`` (``PositionSpecificScoringMatrix.mean`` / ``.std``); ``max`` and ``min`` are
    the scores of the consensus and the anticonsensus.
    """

    motif = _motif(instances, counts, pseudocounts, background)
    pssm = _pssm(motif)
    chosen = _base_map(background, "background", positive_sum=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mean = pssm.mean(chosen)
        std = pssm.std(chosen)
        high, low = pssm.max, pssm.min
    return {
        "status": "completed",
        "operation": "pssm_score_statistics",
        "width": int(motif.length),
        "mean": _finite(mean),
        "std": _finite(std),
        "max": _finite(high),
        "min": _finite(low),
        "provenance": _provenance(
            parser="Bio.motifs.matrix.PositionSpecificScoringMatrix.mean / std / max / min"
        ),
    }


def _distribution(motif: Any, precision: int, background: Any) -> Any:
    width = int(motif.length)
    if precision * width > _MAX_DISTRIBUTION_POINTS:
        raise ConfigurationError(
            f"precision * motif width must be at most {_MAX_DISTRIBUTION_POINTS}; lower the precision."
        )
    pssm = _pssm(motif)
    chosen = _base_map(background, "background", positive_sum=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return pssm.distribution(background=chosen, precision=precision)


def pssm_score_distribution(
    instances: Sequence[str] | None = None,
    *,
    counts: Mapping[str, Sequence[float]] | None = None,
    pseudocounts: float | Mapping[str, float] | None = 0.5,
    background: Mapping[str, float] | None = None,
    precision: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Approximate score distribution of a PSSM for the motif and for the background.

    This is ``PositionSpecificScoringMatrix.distribution`` (``Bio.motifs.thresholds.ScoreDistribution``):
    a dynamic programme over ``precision * width`` score bins. Returns the score of each bin
    (``min_score + i * step``) with the probability of that score under the motif model
    (``motif_density``) and under the background (``background_density``).
    """

    motif = _motif(instances, counts, pseudocounts, background)
    count = _integer(precision, "precision", minimum=1, maximum=100_000)
    dist = _distribution(motif, count, background)
    return {
        "status": "completed",
        "operation": "pssm_score_distribution",
        "precision": count,
        "min_score": float(dist.min_score),
        "step": float(dist.step),
        "bin_count": int(dist.n_points),
        "scores": [float(dist.min_score + index * dist.step) for index in range(dist.n_points)],
        "motif_density": [float(item) for item in dist.mo_density],
        "background_density": [float(item) for item in dist.bg_density],
        "provenance": _provenance(parser="Bio.motifs.thresholds.ScoreDistribution"),
    }


def pssm_threshold(
    instances: Sequence[str] | None = None,
    *,
    counts: Mapping[str, Sequence[float]] | None = None,
    method: str = "fpr",
    rate: float | None = None,
    rate_proportion: float = 1.0,
    pseudocounts: float | Mapping[str, float] | None = 0.5,
    background: Mapping[str, float] | None = None,
    precision: int = 1000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Choose a score threshold for a motif PSSM (``Bio.motifs.thresholds.ScoreDistribution``).

    ``method`` is ``fpr`` (threshold for a false-positive rate ``rate`` against the background),
    ``fnr`` (threshold for a false-negative rate ``rate`` against the motif), ``balanced`` (the threshold
    where the false-negative rate is ``rate_proportion`` times the false-positive rate; also returns the
    rate) or ``patser`` (the threshold used by the Patser program, a balanced threshold from the
    information content). ``rate`` is required for ``fpr`` and ``fnr``.
    """

    chosen = method if method in _THRESHOLD_METHODS else None
    if chosen is None:
        raise ConfigurationError(f"method must be one of {list(_THRESHOLD_METHODS)}; got {method!r}.")
    motif = _motif(instances, counts, pseudocounts, background)
    dist = _distribution(motif, _integer(precision, "precision", minimum=1, maximum=100_000), background)
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "pssm_threshold",
        "method": chosen,
        "precision": precision,
    }
    if chosen in ("fpr", "fnr"):
        if isinstance(rate, bool) or not isinstance(rate, (int, float)) or not 0 < rate < 1:
            raise ConfigurationError("rate must be a number between 0 and 1 for methods fpr and fnr.")
        result["rate"] = float(rate)
        result["threshold"] = float(dist.threshold_fpr(rate) if chosen == "fpr" else dist.threshold_fnr(rate))
    elif chosen == "balanced":
        if isinstance(rate_proportion, bool) or not isinstance(rate_proportion, (int, float)) or rate_proportion <= 0:
            raise ConfigurationError("rate_proportion must be a positive number.")
        threshold, achieved = dist.threshold_balanced(rate_proportion=rate_proportion, return_rate=True)
        result.update({"rate_proportion": float(rate_proportion), "threshold": float(threshold), "rate": float(achieved)})
    else:
        result["threshold"] = float(dist.threshold_patser())
    result["provenance"] = _provenance(
        parser={
            "fpr": "Bio.motifs.thresholds.ScoreDistribution.threshold_fpr",
            "fnr": "Bio.motifs.thresholds.ScoreDistribution.threshold_fnr",
            "balanced": "Bio.motifs.thresholds.ScoreDistribution.threshold_balanced",
            "patser": "Bio.motifs.thresholds.ScoreDistribution.threshold_patser",
        }[chosen]
    )
    return result


def jaspar_pseudocounts(
    instances: Sequence[str] | None = None,
    *,
    counts: Mapping[str, Sequence[float]] | None = None,
    background: Mapping[str, float] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pseudocounts as JASPAR computes them: square root of the number of sites times the background.

    This is ``Bio.motifs.jaspar.calculate_pseudocounts``. The number of sites is the column total of the
    count matrix (the largest column total if columns differ). Returns one pseudocount per base.
    """

    motif = _motif(instances, counts, None, background)
    jaspar = _bio("Bio.motifs.jaspar")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        values = jaspar.calculate_pseudocounts(motif)
    return {
        "status": "completed",
        "operation": "jaspar_pseudocounts",
        "pseudocounts": {base: float(values[base]) for base in _BASES},
        "provenance": _provenance(parser="Bio.motifs.jaspar.calculate_pseudocounts"),
    }


def weblogo_service_logo(
    instances: Sequence[str] | None,
    output_path: str | Path,
    *,
    counts: Mapping[str, Sequence[float]] | None = None,
    image_format: str = "png",
    options: Mapping[str, str | bool] | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw a sequence logo with the WebLogo web service (``Bio.motifs.Motif.weblogo``).

    The motif is sent to https://weblogo.threeplusone.com (WebLogo 3.3) and the returned image is saved
    to a new file ``output_path``. ``image_format`` is png, png_print, pdf, jpeg, eps, svg or txt.
    ``options`` are WebLogo form fields (for example ``{"logo_title": "LexA", "unit_name": "probability",
    "show_errorbars": False}``); unknown field names are refused. Needs network access. Local logos
    without a web service are made by ``plot_pwm_logo``.
    """

    selected = image_format.lower() if isinstance(image_format, str) else ""
    if selected not in _LOGO_FORMATS:
        raise ConfigurationError(f"image_format must be one of {list(_LOGO_FORMATS)}.")
    extra: dict[str, Any] = {}
    for key, value in (options or {}).items():
        expected = _LOGO_OPTIONS.get(key)
        if expected is None:
            raise ConfigurationError(f"{key!r} is not a WebLogo option; allowed: {sorted(_LOGO_OPTIONS)}.")
        if not isinstance(value, expected) or (expected is str and (len(value) > 200 or "\n" in value)):
            raise ConfigurationError(f"option {key} must be {expected.__name__}.")
        extra[key] = value
    motif = _motif(instances, counts, None, None)
    target = cli.new_file(output_path, "output_path")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            motif.weblogo(str(target), fmt=selected.upper(), **extra)
    except Exception as exc:  # urllib raises several network error types
        target.unlink(missing_ok=True)
        raise BackendExecutionError(f"The WebLogo service request failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "weblogo_service_logo",
        "service": "https://weblogo.threeplusone.com/create.cgi (WebLogo 3.3)",
        "output": str(target),
        "image_format": selected,
        "bytes": target.stat().st_size,
        "width": int(motif.length),
        "provenance": _provenance(parser="Bio.motifs.Motif.weblogo"),
    }


__all__ = [
    "jaspar_pseudocounts",
    "motif_anticonsensus",
    "motif_matrices",
    "pssm_score_distribution",
    "pssm_score_statistics",
    "pssm_search",
    "pssm_threshold",
    "weblogo_service_logo",
]

_register_public_verbose(globals(), tuple(__all__))
