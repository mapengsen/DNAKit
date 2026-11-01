"""Motif comparison, matrix transformation, export and scan functions over gimmemotifs (pack P28).

Every function builds ``gimmemotifs.motif.Motif`` objects from explicit PFM count matrices (one
row per position, columns A, C, G and T) or from a PFM file and calls the gimmemotifs routine of
the same name. Results are plain JSON-serialisable dictionaries. Nothing is downloaded.
"""

from __future__ import annotations

import math
import tempfile
from collections.abc import Sequence
from pathlib import Path
from types import SimpleNamespace
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, require
from dnakit.exceptions import BackendExecutionError, ConfigurationError
from dnakit.regulation.backends import (
    _choose,
    _flag,
    _gimme_provenance,
    _integer,
    _pfm,
)

_GIMMEMOTIFS_EXTRA = "bio-gimmemotifs"
_COMPARE_METRICS = ("wic", "pcc", "ed", "distance", "chisq", "akl", "ssd", "seqcor")
_PVALUE_METRICS = ("wic", "pcc", "ed", "distance")
_MATCHES = ("total", "subtotal", "partial")
_COMBINE = ("mean", "sum")
_FORMATS = ("ppm", "pfm", "meme", "motevo", "transfac")
_SCAN_OUTPUTS = ("scores", "positions", "hits", "gff")
_MAX_SEQUENCES = 100_000
_MAX_SAMPLES = 1_000_000
_MAX_MOTIF_FILE_ROWS = 1_000


def _call(function: Any, *args: Any, **kwargs: Any) -> Any:
    """Run an upstream call and convert its failures into DNAKit errors."""

    try:
        return function(*args, **kwargs)
    except KeyboardInterrupt:
        raise
    except BaseException as exc:  # gimmemotifs raises several non-DNAKit exception types
        raise BackendExecutionError(f"gimmemotifs failed: {type(exc).__name__}: {exc}") from exc


def _motif_class() -> Any:
    return require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA).Motif


def _motif(pfm: object, name: str, motif_id: str | None = None) -> Any:
    rows = [[int(v) if v.is_integer() else v for v in row] for row in _pfm(pfm, name)]
    motif = _call(_motif_class(), pfm=rows)
    if motif_id is not None:
        motif.id = motif_id
    return motif


def _finite(value: object) -> float | None:
    if value is None:
        return None
    number = float(value)  # type: ignore[arg-type]
    return number if math.isfinite(number) else None


def _matrix(array: Any, digits: int | None = None) -> list[list[float]]:
    rows = [[float(v) for v in row] for row in array]
    if digits is not None:
        rows = [[round(v, digits) for v in row] for row in rows]
    return rows


def _seed(seed: object) -> int | None:
    if seed is None:
        return None
    return _integer(seed, "seed", 0, 2**32 - 1)


def gimmemotifs_compare_motif_pair(
    pfm_a: Sequence[Sequence[float]],
    pfm_b: Sequence[Sequence[float]],
    *,
    metric: str = "wic",
    match: str = "total",
    combine: str = "mean",
    pvalue: bool = False,
    score_distribution: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Compare two motifs over all alignments with ``MotifComparer.compare_motifs``.

    ``metric`` is wic (weighted information content, the gimmemotifs default), pcc, ed, distance,
    chisq, akl, ssd or seqcor. ``match`` is ``total`` (the whole match, including unaligned
    positions), ``subtotal`` or ``partial`` (only the overlapping positions); ``combine`` is
    ``mean`` or ``sum``. With ``pvalue=True`` the score is converted to the p-value of a normal
    distribution fitted per motif-length pair (``MotifComparer.pvalue``; match total or partial).
    gimmemotifs bundles that distribution for ``metric="wic"``, ``match="total"`` and
    ``combine="mean"`` only; for any other metric pass ``score_distribution``, a file written by
    ``gimmemotifs_motif_score_distribution`` (match ``total``). The result holds the best score
    (or p-value), the offset of the second motif against the first and its strand (1 or -1).
    """

    chosen_metric = _choose(metric, _COMPARE_METRICS, "metric")
    chosen_match = _choose(match, _MATCHES, "match")
    chosen_combine = _choose(combine, _COMBINE, "combine")
    want_p = _flag(pvalue, "pvalue")
    if want_p and chosen_metric not in _PVALUE_METRICS:
        raise ConfigurationError(f"pvalue is available for {', '.join(_PVALUE_METRICS)} only.")
    if score_distribution is not None and not want_p:
        raise ConfigurationError("score_distribution is only used with pvalue=True.")
    if want_p and chosen_match == "subtotal":
        raise ConfigurationError("pvalue is available for match total or partial only.")
    comparison = require("gimmemotifs.comparison", _GIMMEMOTIFS_EXTRA)
    first = _motif(pfm_a, "pfm_a", "motif_a")
    second = _motif(pfm_b, "pfm_b", "motif_b")
    comparer = _call(comparison.MotifComparer)
    if want_p:
        if score_distribution is not None:
            table: dict[int, dict[int, list[float]]] = {}
            for line in input_path(score_distribution).read_text().splitlines():
                length_a, length_b, mean, sd = line.strip().split("\t")[:4]
                table.setdefault(int(length_a), {})[int(length_b)] = [float(mean), float(sd)]
            comparer.scoredist[chosen_metric][f"total_{chosen_combine}"] = table
        if not comparer.scoredist.get(chosen_metric, {}).get(f"total_{chosen_combine}"):
            raise ConfigurationError(
                f"No score distribution for metric {chosen_metric} and combine {chosen_combine}; "
                "create one with gimmemotifs_motif_score_distribution and pass it as "
                "score_distribution."
            )
    outcome = _call(
        comparer.compare_motifs,
        first,
        second,
        chosen_match,
        chosen_metric,
        chosen_combine,
        pval=want_p,
    )
    if outcome is None or len(outcome) != 3:
        raise BackendExecutionError("gimmemotifs returned no comparison result.")
    value, position, strand = outcome
    return {
        "status": "completed",
        "operation": "compare_motif_pair",
        "metric": chosen_metric,
        "match": chosen_match,
        "combine": chosen_combine,
        "value_kind": "pvalue" if want_p else "score",
        "value": _finite(value),
        "position": None
        if position is None or not math.isfinite(float(position))
        else int(position),
        "strand": None if strand is None or not math.isfinite(float(strand)) else int(strand),
        "length_a": len(first),
        "length_b": len(second),
        "provenance": _gimme_provenance(
            function="gimmemotifs.comparison.MotifComparer.compare_motifs"
        ),
    }


def gimmemotifs_motif_score_distribution(
    motifs_pfm: str | Path,
    output_file: str | Path,
    *,
    metric: str = "wic",
    match: str = "total",
    combine: str = "mean",
    verbose: bool = True,
) -> dict[str, Any]:
    """Score distribution of all motif pairs by length, as used for motif-similarity p-values.

    Runs ``MotifComparer.generate_score_dist`` on the motifs of a PFM file: every motif is compared
    with every motif and, for each pair of motif lengths, the mean and standard deviation of the
    scores are written as ``length1<TAB>length2<TAB>mean<TAB>sd`` lines to ``output_file``. Upstream
    writes into the gimmemotifs installation; here the file goes to ``output_file`` (must not
    exist) and the installed score tables stay untouched. ``metric`` is wic, pcc, ed or distance.
    """

    chosen_metric = _choose(metric, _PVALUE_METRICS, "metric")
    chosen_match = _choose(match, ("total", "subtotal"), "match")
    chosen_combine = _choose(combine, _COMBINE, "combine")
    source = input_path(motifs_pfm)
    destination = Path(output_file).expanduser()
    if destination.exists():
        raise ConfigurationError(f"output_file already exists: {destination}")
    comparison = require("gimmemotifs.comparison", _GIMMEMOTIFS_EXTRA)
    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    motifs = _call(motif_module.read_motifs, str(source))
    if not 1 <= len(motifs) <= _MAX_MOTIF_FILE_ROWS:
        raise ConfigurationError(f"motifs_pfm must hold 1 to {_MAX_MOTIF_FILE_ROWS} motifs.")
    comparer = _call(comparison.MotifComparer)
    original = comparer.get_all_scores

    def serial(*args: Any, **kwargs: Any) -> Any:
        kwargs["parallel"] = False
        return original(*args, **kwargs)

    with tempfile.TemporaryDirectory(prefix="dnakit_gimme_scoredist_") as workspace:
        comparer.config = SimpleNamespace(get_score_dir=lambda: workspace)
        comparer.get_all_scores = serial
        _call(comparer.generate_score_dist, motifs, chosen_match, chosen_metric, chosen_combine)
        written = (
            Path(workspace) / f"{chosen_match}_{chosen_metric}_{chosen_combine}_score_dist.txt"
        )
        text = written.read_text()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(text)
    rows = []
    for line in text.splitlines():
        first, second, mean, sd = line.split("\t")
        rows.append(
            {"length_a": int(first), "length_b": int(second), "mean": float(mean), "sd": float(sd)}
        )
    return {
        "status": "completed",
        "operation": "motif_score_distribution",
        "metric": chosen_metric,
        "match": chosen_match,
        "combine": chosen_combine,
        "n_motifs": len(motifs),
        "distribution": rows,
        "outputs": {"score_distribution": str(destination)},
        "provenance": _gimme_provenance(
            function="gimmemotifs.comparison.MotifComparer.generate_score_dist"
        ),
    }


def gimmemotifs_motif_pcc(
    pfm_a: Sequence[Sequence[float]],
    pfm_b: Sequence[Sequence[float]],
    *,
    position: int | None = None,
    revcomp: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Pearson correlation of two motifs' position probabilities.

    Calls ``Motif.pcc`` and ``Motif.max_pcc``.

    With ``position`` (an integer offset of the second motif against the first; -2 puts the second
    motif two positions to the left) the sum of the per-position Pearson correlations of the
    overlapping positions is returned (``pcc``). Without ``position`` every alignment is tried,
    also with the reverse complement of the second motif when ``revcomp`` is true, and the best
    one is returned (``max_pcc``) with its offset and strand (1 or -1).
    """

    first = _motif(pfm_a, "pfm_a", "motif_a")
    second = _motif(pfm_b, "pfm_b", "motif_b")
    _flag(revcomp, "revcomp")
    if position is not None:
        offset = _integer(position, "position", -10_000, 10_000)
        score = _call(first.pcc, first.ppm, second.ppm, offset)
        return {
            "status": "completed",
            "operation": "motif_pcc",
            "mode": "pcc",
            "score": _finite(score),
            "position": offset,
            "strand": 1,
            "provenance": _gimme_provenance(function="gimmemotifs.motif._comparison.pcc"),
        }
    score, offset, strand = _call(first.max_pcc, second, revcomp=revcomp)
    return {
        "status": "completed",
        "operation": "motif_pcc",
        "mode": "max_pcc",
        "score": _finite(score),
        "position": int(offset),
        "strand": int(strand),
        "provenance": _gimme_provenance(function="gimmemotifs.motif._comparison.max_pcc"),
    }


def gimmemotifs_motif_ic_similarity(
    pfm_a: Sequence[Sequence[float]],
    pfm_b: Sequence[Sequence[float]],
    *,
    mode: str = "matrix",
    position: int = 0,
    bg_factor: float = 1.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Information-content similarity of two motifs (``Motif.matrix_ic`` and ``Motif.ic``).

    The result also carries the information content of each motif (``Motif.information_content``,
    bits). ``mode="matrix"`` tries every alignment of both strands and returns the best score, its
    offset and strand (``matrix_ic``). ``mode="aligned"`` scores the single alignment at
    ``position`` (``ic``; unaligned positions are compared with the uniform background, weighted
    by ``bg_factor``). GimmeMotifs ``matrix_ic`` raises a broadcasting error for some pairs in
    which the second motif is longer than the first (an upstream defect); the error is reported
    with a hint.
    """

    chosen = _choose(mode, ("matrix", "aligned"), "mode")
    first = _motif(pfm_a, "pfm_a", "motif_a")
    second = _motif(pfm_b, "pfm_b", "motif_b")
    offset = _integer(position, "position", -10_000, 10_000)
    factor = float(bg_factor)
    if not math.isfinite(factor) or factor < 0:
        raise ConfigurationError("bg_factor must be a finite non-negative number.")
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "motif_ic_similarity",
        "mode": chosen,
        "information_content_a": _finite(first.information_content),
        "information_content_b": _finite(second.information_content),
    }
    if chosen == "matrix":
        try:
            score, best, strand = _call(first.matrix_ic, first.ppm, second.ppm)
        except BackendExecutionError as exc:
            raise BackendExecutionError(
                f"{exc} Motif.matrix_ic fails for some pairs of motif lengths (upstream defect "
                "when the second motif is longer than the first); swap the motifs or use "
                "mode='aligned'."
            ) from exc
        result.update(score=_finite(score), position=int(best), strand=int(strand))
        result["provenance"] = _gimme_provenance(function="gimmemotifs.motif._comparison.matrix_ic")
    else:
        score = _call(first.ic, first.ppm.tolist(), second.ppm.tolist(), offset, bg_factor=factor)
        result.update(score=_finite(score), position=offset, strand=1)
        result["provenance"] = _gimme_provenance(function="gimmemotifs.motif._comparison.ic")
    return result


def gimmemotifs_motif_score_range(
    pfm: Sequence[Sequence[float]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Minimum and maximum PWM.

    Calls log-odds) score of a motif (``Motif.pwm_min_score`` and ``pwm_max_score``.

    The scores are the sums over positions of the smallest and largest log-odds value
    ``log(p / 0.25 + 0.01)`` that gimmemotifs scans with. Also returns the motif length.
    """

    motif = _motif(pfm, "pfm", "motif")
    return {
        "status": "completed",
        "operation": "motif_score_range",
        "length": len(motif),
        "min_score": _finite(motif.min_score),
        "max_score": _finite(motif.max_score),
        "provenance": _gimme_provenance(
            function="gimmemotifs.motif.base.Motif.min_score/max_score"
        ),
    }


def gimmemotifs_motif_consensus(
    pfm: Sequence[Sequence[float]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """IUPAC consensus sequence of a motif (``Motif.to_consensusv2``).

    Upper-case letters mark positions dominated by one or two bases, lower-case letters weaker
    ambiguity codes and ``n`` an uninformative position.
    """

    motif = _motif(pfm, "pfm", "motif")
    consensus = _call(motif.to_consensusv2)
    return {
        "status": "completed",
        "operation": "motif_consensus",
        "consensus": str(consensus),
        "length": len(motif),
        "provenance": _gimme_provenance(function="gimmemotifs.motif.base.Motif.to_consensusv2"),
    }


def gimmemotifs_motif_reverse_complement(
    pfm: Sequence[Sequence[float]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Reverse-complement a motif matrix (``Motif.rc``): reverse positions and A/C/G/T columns."""

    motif = _motif(pfm, "pfm", "motif")
    reverse = _call(motif.rc)
    return {
        "status": "completed",
        "operation": "motif_reverse_complement",
        "pfm": _matrix(reverse.pfm),
        "ppm": _matrix(reverse.ppm),
        "consensus": str(reverse.consensus),
        "provenance": _gimme_provenance(function="gimmemotifs.motif.base.Motif.rc"),
    }


def gimmemotifs_motif_shuffle(
    pfm: Sequence[Sequence[float]],
    *,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Shuffle the positions of a motif (``Motif.shuffle``), a standard control motif.

    Upstream draws from numpy's global random state; ``seed`` seeds it before the call so that the
    permutation is reproducible. The order of the A/C/G/T columns is kept.
    """

    np = require("numpy", _GIMMEMOTIFS_EXTRA)
    motif = _motif(pfm, "pfm", "motif")
    chosen = _seed(seed)
    state = np.random.get_state()
    try:
        if chosen is not None:
            np.random.seed(chosen)
        shuffled = _call(motif.shuffle)
    finally:
        np.random.set_state(state)
    return {
        "status": "completed",
        "operation": "motif_shuffle",
        "seed": chosen,
        "pfm": _matrix(shuffled.pfm),
        "ppm": _matrix(shuffled.ppm),
        "provenance": _gimme_provenance(function="gimmemotifs.motif.base.Motif.shuffle"),
    }


def gimmemotifs_pfm_to_ppm(
    pfm: Sequence[Sequence[float]],
    *,
    pseudo: float = 0.001,
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert motif counts to frequencies with a pseudocount (``Motif.pfm_to_ppm``).

    Every cell becomes ``(count + pseudo) / (row sum + 4 * pseudo)``, so each row sums to 1.
    """

    rows = _pfm(pfm, "pfm")
    value = float(pseudo)
    if not math.isfinite(value) or value < 0:
        raise ConfigurationError("pseudo must be a finite non-negative number.")
    motif = _call(_motif_class())
    ppm = _call(motif.pfm_to_ppm, rows, pseudo=value)
    return {
        "status": "completed",
        "operation": "pfm_to_ppm",
        "pseudo": value,
        "ppm": _matrix(ppm),
        "provenance": _gimme_provenance(function="gimmemotifs.motif.base.Motif.pfm_to_ppm"),
    }


def gimmemotifs_ppm_pseudocount(
    ppm: Sequence[Sequence[float]],
    *,
    pseudo: float = 1e-6,
    verbose: bool = True,
) -> dict[str, Any]:
    """Add a pseudocount to a probability matrix and renormalise each row.

    Calls ``utils.ppm_pseudocount``.

    The default ``1e-6`` is the value gimmemotifs applies before comparing motifs, so no
    probability is exactly zero.
    """

    np = require("numpy", _GIMMEMOTIFS_EXTRA)
    utils = require("gimmemotifs.utils", _GIMMEMOTIFS_EXTRA)
    rows = _pfm(ppm, "ppm")
    value = float(pseudo)
    if not math.isfinite(value) or value < 0:
        raise ConfigurationError("pseudo must be a finite non-negative number.")
    result = _call(utils.ppm_pseudocount, np.asarray(rows, dtype=float), value)
    return {
        "status": "completed",
        "operation": "ppm_pseudocount",
        "pseudo": value,
        "ppm": _matrix(result),
        "provenance": _gimme_provenance(function="gimmemotifs.utils.ppm_pseudocount"),
    }


def gimmemotifs_format_motif(
    pfm: Sequence[Sequence[float]],
    *,
    output_format: str = "ppm",
    motif_id: str = "motif",
    precision: int = 4,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write a motif as text in a gimmemotifs export format.

    ``output_format`` is ``ppm`` (``Motif.to_ppm``: FASTA-like header and tab-separated
    probabilities with ``precision`` decimals), ``pfm`` (``Motif.to_pfm``: counts), ``meme``
    (``Motif.to_meme``), ``motevo`` (``Motif.to_motevo``) or ``transfac`` (``Motif.to_transfac``).
    Returns the text.
    """

    chosen = _choose(output_format, _FORMATS, "output_format")
    name = str(motif_id)
    if not name.strip() or "\n" in name or "\x00" in name:
        raise ConfigurationError("motif_id must be a non-empty single line.")
    digits = _integer(precision, "precision", 0, 12)
    motif = _motif(pfm, "pfm", name)
    if chosen == "ppm":
        text = _call(motif.to_ppm, precision=digits)
    elif chosen == "pfm":
        text = _call(motif.to_pfm)
    elif chosen == "meme":
        text = _call(motif.to_meme)
    elif chosen == "motevo":
        text = _call(motif.to_motevo)
    else:
        text = _call(motif.to_transfac)
    return {
        "status": "completed",
        "operation": "format_motif",
        "output_format": chosen,
        "motif_id": name,
        "text": str(text),
        "provenance": _gimme_provenance(function=f"gimmemotifs.motif.base.Motif.to_{chosen}"),
    }


def gimmemotifs_average_motifs(
    pfm_a: Sequence[Sequence[float]],
    pfm_b: Sequence[Sequence[float]],
    *,
    position: int = 0,
    orientation: int = 1,
    include_bg: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Average two motifs at a given alignment (``Motif.average_motifs``).

    ``position`` is the offset of the second motif against the first (-2 places the second two
    positions to the left). ``orientation`` 1 uses the second motif as given, -1 its reverse
    complement. Positions covered by one motif only are filled with the other motif's
    frequencies, or with the uniform background when ``include_bg`` is true. The counts of the
    two matrices are added.
    """

    first = _motif(pfm_a, "pfm_a", "motif_a")
    second = _motif(pfm_b, "pfm_b", "motif_b")
    offset = _integer(position, "position", -10_000, 10_000)
    if orientation not in (1, -1) or isinstance(orientation, bool):
        raise ConfigurationError("orientation must be 1 or -1.")
    averaged = _call(
        first.average_motifs,
        second,
        offset,
        orientation,
        include_bg=_flag(include_bg, "include_bg"),
    )
    return {
        "status": "completed",
        "operation": "average_motifs",
        "position": offset,
        "orientation": orientation,
        "pfm": _matrix(averaged.pfm),
        "ppm": _matrix(averaged.ppm),
        "consensus": str(averaged.consensus),
        "provenance": _gimme_provenance(function="gimmemotifs.motif.base.Motif.average_motifs"),
    }


def gimmemotifs_sample_motif_sequences(
    pfm: Sequence[Sequence[float]],
    n_sequences: int,
    *,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Sample random sequences that follow a motif's probability matrix (``Motif.sample``).

    Each position is drawn independently from the motif's probabilities. ``seed`` makes the draw
    reproducible (numpy ``default_rng``).
    """

    np = require("numpy", _GIMMEMOTIFS_EXTRA)
    count = _integer(n_sequences, "n_sequences", 1, _MAX_SAMPLES)
    motif = _motif(pfm, "pfm", "motif")
    chosen = _seed(seed)
    rng = np.random.default_rng(chosen)
    sequences = _call(motif.sample, count, rng=rng)
    return {
        "status": "completed",
        "operation": "sample_motif_sequences",
        "seed": chosen,
        "n_sequences": len(sequences),
        "sequences": [str(s) for s in sequences],
        "provenance": _gimme_provenance(function="gimmemotifs.motif.base.Motif.sample"),
    }


def _scan_fasta(fasta: str | Path) -> Any:
    fasta_module = require("gimmemotifs.fasta", _GIMMEMOTIFS_EXTRA)
    path = input_path(fasta)
    fa = _call(fasta_module.Fasta, str(path))
    if not 1 <= len(fa) <= _MAX_SEQUENCES:
        raise ConfigurationError(f"FASTA must hold 1 to {_MAX_SEQUENCES} sequences.")
    return fa


def gimmemotifs_pwm_scan(
    sequences_fasta: str | Path,
    pfm: Sequence[Sequence[float]],
    *,
    output_format: str = "scores",
    cutoff: float | None = None,
    nreport: int | None = None,
    scan_rc: bool = True,
    gff_path: str | Path | None = None,
    motif_id: str = "motif",
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan FASTA sequences with one motif matrix.

    Calls ``Motif.pwm_scan_score``, ``scan_all``, ``scan``, ``pwm_scan_to_gff``.

    ``cutoff`` is the fraction (0 to 1) of the way from the minimum to the maximum PWM score
    that a hit must reach. ``output_format`` selects the output: ``scores`` (the score of each
    reported match per sequence; defaults ``cutoff=0``, ``nreport=1`` = best match, as upstream),
    ``hits`` (position, score and strand per match; defaults ``cutoff=0.9``, ``nreport=50``),
    ``positions`` (match positions only, same defaults) or ``gff`` (also writes a GFF file to
    ``gff_path``, which must not exist; same defaults as ``hits``). ``scan_rc`` includes the
    reverse strand. Sequence names are the complete FASTA header lines.
    """

    chosen = _choose(output_format, _SCAN_OUTPUTS, "output_format")
    default_cutoff, default_nreport = (0.0, 1) if chosen == "scores" else (0.9, 50)
    fraction = default_cutoff if cutoff is None else float(cutoff)
    if not math.isfinite(fraction) or not 0.0 <= fraction <= 1.0:
        raise ConfigurationError("cutoff must be a number in [0, 1].")
    limit = default_nreport if nreport is None else _integer(nreport, "nreport", 1, 100_000)
    strand = _flag(scan_rc, "scan_rc")
    motif = _motif(pfm, "pfm", str(motif_id))
    fa = _scan_fasta(sequences_fasta)
    outputs: dict[str, str] = {}
    if chosen == "scores":
        raw = _call(motif.pwm_scan_score, fa, cutoff=fraction, nreport=limit, scan_rc=strand)
        matches: dict[str, Any] = {k: [float(v) for v in vals] for k, vals in raw.items()}
        function = "pwm_scan_score"
    elif chosen == "positions":
        raw = _call(motif.scan, fa, cutoff=fraction, nreport=limit, scan_rc=strand)
        matches = {k: [int(v) for v in vals] for k, vals in raw.items()}
        function = "scan"
    else:
        raw = _call(motif.scan_all, fa, cutoff=fraction, nreport=limit, scan_rc=strand)
        matches = {
            k: [{"position": int(p), "score": float(s), "strand": int(t)} for p, s, t in vals]
            for k, vals in raw.items()
        }
        function = "scan_all"
        if chosen == "gff":
            if gff_path is None:
                raise ConfigurationError("gff_path is required for output_format='gff'.")
            destination = Path(gff_path).expanduser()
            if destination.exists():
                raise ConfigurationError(f"gff_path already exists: {destination}")
            destination.parent.mkdir(parents=True, exist_ok=True)
            _call(
                motif.pwm_scan_to_gff,
                fa,
                str(destination),
                cutoff=fraction,
                nreport=limit,
                scan_rc=strand,
            )
            outputs["gff"] = str(destination)
            function = "pwm_scan_to_gff"
    return {
        "status": "completed",
        "operation": "pwm_scan",
        "output_format": chosen,
        "cutoff": fraction,
        "nreport": limit,
        "scan_rc": strand,
        "n_sequences": len(matches),
        "matches": matches,
        "outputs": outputs,
        "provenance": _gimme_provenance(function=f"gimmemotifs.motif._scanning.{function}"),
    }


def gimmemotifs_consensus_scan(
    sequences_fasta: str | Path,
    pfm: Sequence[Sequence[float]],
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan FASTA sequences for the motif's IUPAC consensus as a regular expression.

    Calls ``Motif.consensus_scan``.

    Each IUPAC symbol of the consensus is expanded into a character class, so only exact
    consensus matches are found on the forward strand. Returns, per sequence, the middle
    coordinate of every match (0-based, half-integers possible).
    """

    motif = _motif(pfm, "pfm", "motif")
    fa = _scan_fasta(sequences_fasta)
    raw = _call(motif.consensus_scan, fa)
    return {
        "status": "completed",
        "operation": "consensus_scan",
        "consensus": str(motif.to_consensusv2()),
        "n_sequences": len(raw),
        "matches": {k: [float(v) for v in vals] for k, vals in raw.items()},
        "provenance": _gimme_provenance(function="gimmemotifs.motif.base.Motif.consensus_scan"),
    }


__all__ = [
    "gimmemotifs_average_motifs",
    "gimmemotifs_compare_motif_pair",
    "gimmemotifs_consensus_scan",
    "gimmemotifs_format_motif",
    "gimmemotifs_motif_consensus",
    "gimmemotifs_motif_ic_similarity",
    "gimmemotifs_motif_pcc",
    "gimmemotifs_motif_reverse_complement",
    "gimmemotifs_motif_score_distribution",
    "gimmemotifs_motif_score_range",
    "gimmemotifs_motif_shuffle",
    "gimmemotifs_pfm_to_ppm",
    "gimmemotifs_ppm_pseudocount",
    "gimmemotifs_pwm_scan",
    "gimmemotifs_sample_motif_sequences",
]

_register_public_verbose(globals(), tuple(__all__))
