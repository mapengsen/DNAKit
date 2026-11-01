"""DnaChisel design rules and PSSM patterns that DNAKit did not yet expose (integration pack P38).

``evaluate_dnachisel_design_rule`` scores (and optionally optimises) a DNA sequence against one DnaChisel design
rule chosen by name: AllowPrimer, AvoidBlastMatches, AvoidHeterodimerization, AvoidMatches, EnforceChoice,
EnforcePatternOccurence, EnforceRegionsCompatibility, EnforceSequence, EnforceTerminalGCContent,
SequenceLengthBounds, UniquifyAllKmers, AvoidRareCodons, MatchTargetCodonUsage and HarmonizeRCA. The two PSSM
functions build DnaChisel ``MotifPssmPattern`` objects from same-length sequences or from motif files.

Install with ``python -m pip install "dnakit[bio-dnachisel]"``. AvoidBlastMatches needs BLAST+ (``blastn``) and
AvoidMatches needs Bowtie 1 on PATH; a missing program raises ``BackendUnavailableError`` from DnaChisel's own
check. Coordinates are 0-based, half-open. Tested with dnachisel 3.2.16.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p38_helpers as h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendUnavailableError, ConfigurationError

_EXTRA = "bio-dnachisel"
_MAX_DNA = 1_000_000
_RULES: dict[str, tuple[str, ...]] = {
    "AllowPrimer": ("location", "tmin", "tmax", "max_homology_length", "avoid_heterodim_with",
                    "max_heterodim_tm", "avoided_repeats", "boost"),
    "AvoidBlastMatches": ("blast_db", "sequences", "word_size", "perc_identity", "num_alignments",
                          "num_threads", "min_align_length", "ungapped", "e_value", "culling_limit",
                          "location", "boost"),
    "AvoidHeterodimerization": ("other_primers_sequences", "tmax", "location", "boost"),
    "AvoidMatches": ("match_length", "bowtie_index", "sequences", "mismatches", "location", "boost"),
    "EnforceChoice": ("choices", "location", "boost"),
    "EnforcePatternOccurence": ("pattern", "occurences", "location", "strand", "center", "boost"),
    "EnforceRegionsCompatibility": ("locations", "condition", "min_distance", "boost"),
    "EnforceSequence": ("sequence", "location", "boost"),
    "EnforceTerminalGCContent": ("window_size", "mini", "maxi", "boost"),
    "SequenceLengthBounds": ("min_length", "max_length", "boost"),
    "UniquifyAllKmers": ("k", "reference", "location", "include_reverse_complement", "boost"),
    "AvoidRareCodons": ("min_frequency", "species", "codon_usage_table", "location", "boost"),
    "MatchTargetCodonUsage": ("species", "location", "codon_usage_table", "boost"),
    "HarmonizeRCA": ("species", "codon_usage_table", "original_species", "original_codon_usage_table",
                     "location", "boost"),
}
_CONSTRAINT_RULES = frozenset({
    "EnforceChoice", "EnforceSequence", "SequenceLengthBounds", "EnforceTerminalGCContent",
    "EnforcePatternOccurence", "EnforceRegionsCompatibility", "AvoidBlastMatches", "AvoidMatches",
    "AllowPrimer", "AvoidHeterodimerization",
})
_CONDITIONS = ("different", "min_hamming_distance")


def _dna(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    text = value.strip().upper()
    if len(text) > _MAX_DNA or set(text) - set("ACGTN"):
        raise ConfigurationError(f"{name} must contain only A, C, G, T or N and at most {_MAX_DNA} bases.")
    return text


def _location(dnachisel: Any, value: Any, name: str = "location") -> Any:
    if value is None:
        return None
    if not isinstance(value, (list, tuple)) or len(value) not in (2, 3):
        raise ConfigurationError(f"{name} must be [start, end] or [start, end, strand].")
    return dnachisel.Location.from_data(tuple(value))


def _build_rule(dnachisel: Any, rule: str, params: Mapping[str, Any]) -> Any:
    allowed = _RULES[rule]
    unknown = sorted(set(params) - set(allowed))
    if unknown:
        raise ConfigurationError(f"Unknown parameter(s) for {rule}: {', '.join(unknown)}. Allowed: {', '.join(allowed)}.")
    kwargs = dict(params)
    if "location" in kwargs:
        kwargs["location"] = _location(dnachisel, kwargs["location"])
    if rule == "AllowPrimer" and "avoided_repeats" in kwargs:
        kwargs["avoided_repeats"] = tuple(tuple(int(v) for v in pair) for pair in kwargs["avoided_repeats"])
    if rule == "EnforceRegionsCompatibility":
        locations = kwargs.pop("locations", None)
        if not isinstance(locations, (list, tuple)) or len(locations) < 2:
            raise ConfigurationError("locations must list at least two [start, end] regions.")
        condition = h.choice(kwargs.pop("condition", "different"), "condition", _CONDITIONS)
        min_distance = h.integer(kwargs.pop("min_distance", 1), "min_distance", 1, 10_000)
        if condition == "different":
            def compatible(first: Any, second: Any, problem: Any) -> bool:
                return bool(first.extract_sequence(problem.sequence) != second.extract_sequence(problem.sequence))
            label = "different sequences "
        else:
            def compatible(first: Any, second: Any, problem: Any) -> bool:
                a, b = first.extract_sequence(problem.sequence), second.extract_sequence(problem.sequence)
                return len(a) == len(b) and sum(x != y for x, y in zip(a, b)) >= min_distance
            label = f"Hamming distance >= {min_distance} "
        pairs = []
        for loc in locations:
            if not isinstance(loc, (list, tuple)) or len(loc) not in (2, 3):
                raise ConfigurationError("Each locations entry must be [start, end] or [start, end, strand].")
            pairs.append(tuple(int(v) for v in loc))
        return dnachisel.EnforceRegionsCompatibility(
            locations=pairs, compatibility_condition=compatible, condition_label=label, **kwargs,
        )
    try:
        return getattr(dnachisel, rule)(**kwargs)
    except (TypeError, ValueError, KeyError) as exc:
        raise ConfigurationError(f"Invalid parameters for {rule}: {exc}") from exc


def _summary(evaluation: Any) -> dict[str, Any]:
    locations = evaluation.locations
    return {
        "passes": bool(evaluation.passes),
        "score": h.plain(float(evaluation.score)) if evaluation.score is not None else None,
        "message": str(evaluation.message),
        "locations": None if locations is None else [str(loc) for loc in locations],
    }


def evaluate_dnachisel_design_rule(
    sequence: str,
    rule: str,
    *,
    parameters: Mapping[str, Any] | None = None,
    circular: bool = False,
    role: str = "objective",
    optimize: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score a DNA sequence against one DnaChisel design rule, and optionally optimise it for the rule.

    ``rule`` selects the specification: ``AllowPrimer`` (primer-friendly region: Tm window, no repeats, homologies
    or heterodimers), ``AvoidBlastMatches`` (BLAST+ ``blast_db`` or ``sequences``), ``AvoidHeterodimerization``
    (``other_primers_sequences``, ``tmax``), ``AvoidMatches`` (``match_length`` and a Bowtie ``bowtie_index``),
    ``EnforceChoice`` (``choices``), ``EnforcePatternOccurence`` (``pattern``, ``occurences``),
    ``EnforceRegionsCompatibility`` (``locations`` and ``condition`` ``different`` or ``min_hamming_distance`` with
    ``min_distance``), ``EnforceSequence`` (degenerate ``sequence`` at ``location``), ``EnforceTerminalGCContent``
    (``window_size``, ``mini``, ``maxi``), ``SequenceLengthBounds`` (``min_length``, ``max_length``),
    ``UniquifyAllKmers`` (``k``), ``AvoidRareCodons`` (``min_frequency``, ``species``), ``MatchTargetCodonUsage``
    (``species``) and ``HarmonizeRCA`` (``original_species``, ``species``). ``parameters`` holds the DnaChisel
    keyword arguments of that rule; ``location`` is ``[start, end]`` or ``[start, end, strand]``.

    ``role`` is ``objective`` (the sequence is only scored; default) or ``constraint`` (DnaChisel first resolves
    constraints that restrict nucleotides, which can change the sequence). With ``optimize=True`` the rule is used
    as an objective (or constraint) and the sequence is optimised with DnaChisel's random-mutation optimiser; the
    result holds the evaluation before and after and the new sequence. Optimisation is stochastic in DnaChisel.
    """

    name = h.choice(rule, "rule", tuple(_RULES))
    text = _dna(sequence, "sequence")
    if parameters is not None and not isinstance(parameters, Mapping):
        raise ConfigurationError("parameters must be a mapping of keyword arguments.")
    kind = h.choice(role, "role", ("objective", "constraint"))
    use_circular = h.flag(circular, "circular")
    do_optimize = h.flag(optimize, "optimize")
    dnachisel = require("dnachisel", _EXTRA)
    problem_class = dnachisel.CircularDnaOptimizationProblem if use_circular else dnachisel.DnaOptimizationProblem

    def make_problem(seq: str) -> Any:
        spec = _build_rule(dnachisel, name, dict(parameters or {}))
        return problem_class(
            sequence=seq, constraints=[spec] if kind == "constraint" else [],
            objectives=[spec] if kind == "objective" else [], logger=None,
        )

    def evaluate(problem: Any) -> dict[str, Any]:
        items = list(problem.constraints_evaluations() if kind == "constraint" else problem.objectives_evaluations())
        return _summary(items[0])

    try:
        problem = make_problem(text)
        before = evaluate(problem)
        result_sequence = problem.sequence
        after = None
        if do_optimize:
            if kind == "constraint":
                problem.resolve_constraints()
            else:
                problem.optimize()
            result_sequence = problem.sequence
            after = evaluate(problem)
    except BackendUnavailableError:
        raise
    except ConfigurationError:
        raise
    except FileNotFoundError as exc:
        raise BackendUnavailableError(
            f"{name} needs an external program that was not found: {exc}",
            hint="Install BLAST+ (blastn) for AvoidBlastMatches or Bowtie 1 for AvoidMatches and put it on PATH.",
        ) from exc
    except Exception as exc:  # DnaChisel raises NoSolutionError, ValueError, AttributeError, OSError ...
        raise ConfigurationError(f"DnaChisel could not evaluate {name}: {exc}") from exc
    return {
        "status": "completed",
        "operation": "evaluate_dnachisel_design_rule",
        "rule": name,
        "role": kind,
        "circular": use_circular,
        "sequence_length": len(text),
        "evaluation": before if after is None else after,
        "evaluation_before": before,
        "evaluation_after": after,
        "optimized": do_optimize,
        "sequence": result_sequence,
        "sequence_changed": result_sequence != text,
        "constraint_default": name in _CONSTRAINT_RULES,
        "provenance": h.provenance("dnachisel", f"dnachisel.{name}"),
    }


def _pssm_payload(pattern: Any, scan: str | None) -> dict[str, Any]:
    matrix = pattern.pssm
    out: dict[str, Any] = {
        "name": str(pattern.name),
        "length": int(pattern.size),
        "threshold": None if pattern.threshold is None else h.plain(float(pattern.threshold)),
        "relative_threshold": pattern.relative_threshold,
        "min_score": h.plain(float(matrix.min)),
        "max_score": h.plain(float(matrix.max)),
        "pssm": {base: [h.plain(float(v)) for v in matrix[base]] for base in "ACGT"},
    }
    if scan is not None:
        if pattern.threshold is None:
            raise ConfigurationError("Give threshold or relative_threshold to scan a sequence.")
        out["matches"] = [{"start": int(s), "end": int(e), "strand": int(st)}
                          for s, e, st in pattern.find_matches_in_string(scan)]
    return out


def _pseudocounts(value: Any) -> Any:
    if value is None or value == "jaspar":
        return value
    if isinstance(value, Mapping) and set(value) == set("ACGT"):
        return {k: h.number(v, f"pseudocount {k}", 0, 1000) for k, v in value.items()}
    raise ConfigurationError("pseudocounts must be 'jaspar', null or a mapping with A, C, G and T.")


def _thresholds(threshold: float | None, relative_threshold: float | None) -> tuple[float | None, float | None]:
    if threshold is not None and relative_threshold is not None:
        raise ConfigurationError("Give threshold or relative_threshold, not both.")
    t = None if threshold is None else h.number(threshold, "threshold", -1e6, 1e6)
    r = None if relative_threshold is None else h.number(relative_threshold, "relative_threshold", 0, 1)
    return t, r


def build_pssm_pattern_from_sequences(
    sequences: Sequence[str],
    *,
    name: str = "unnamed",
    pseudocounts: Any = "jaspar",
    threshold: float | None = None,
    relative_threshold: float | None = None,
    scan_sequence: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a PSSM pattern from same-length DNA sequences with pseudocounts, using
    ``dnachisel.MotifPssmPattern.from_sequences`` (and ``apply_pseudocounts``).

    ``pseudocounts`` is ``jaspar`` (automatic, Biopython), a mapping ``{"A": .., "C": .., "G": .., "T": ..}`` or
    null for none. A match threshold is given absolute (``threshold``) or relative (``relative_threshold``, 0 to
    1 between the minimum and maximum possible score). With ``scan_sequence`` the forward-strand matches are
    returned. Returns the log-odds matrix by base.
    """

    if isinstance(sequences, (str, bytes)) or not isinstance(sequences, Sequence) or not sequences:
        raise ConfigurationError("sequences must be a non-empty list of equal-length DNA strings.")
    seqs = [_dna(s, "sequence") for s in sequences]
    if len({len(s) for s in seqs}) != 1 or "N" in "".join(seqs):
        raise ConfigurationError("sequences must have the same length and contain only A, C, G and T.")
    t, r = _thresholds(threshold, relative_threshold)
    dnachisel = require("dnachisel", _EXTRA)
    require("Bio.motifs.jaspar", _EXTRA)  # dnachisel's pseudocount routine expects Bio.motifs.jaspar to be loaded
    scan = None if scan_sequence is None else _dna(scan_sequence, "scan_sequence")
    try:
        pattern = dnachisel.MotifPssmPattern.from_sequences(
            seqs, name=h.text(name, "name"), pseudocounts=_pseudocounts(pseudocounts),
            threshold=t, relative_threshold=r,
        )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"DnaChisel could not build the PSSM: {exc}") from exc
    return {
        "status": "completed",
        "operation": "build_pssm_pattern_from_sequences",
        "sequence_count": len(seqs),
        **_pssm_payload(pattern, scan),
        "provenance": h.provenance("dnachisel", "dnachisel.MotifPssmPattern.from_sequences"),
    }


def read_pssm_patterns_from_motif_file(
    motifs_path: str,
    file_format: str,
    *,
    pseudocounts: Any = "jaspar",
    threshold: float | None = None,
    relative_threshold: float | None = None,
    scan_sequence: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read motifs from a JASPAR, MEME or TRANSFAC file and return them as PSSM patterns, using
    ``dnachisel.MotifPssmPattern.list_from_file`` (Bio.motifs parsers).

    ``file_format`` is ``jaspar``, ``meme`` (MEME XML, as read by Biopython) or ``TRANSFAC``. Pseudocounts and thresholds are as in
    ``build_pssm_pattern_from_sequences``; with ``scan_sequence`` the forward-strand matches of every motif are
    listed.
    """

    source = h.existing(motifs_path, "motifs_path")
    fmt = h.choice(file_format, "file_format", ("jaspar", "meme", "TRANSFAC"))
    t, r = _thresholds(threshold, relative_threshold)
    dnachisel = require("dnachisel", _EXTRA)
    require("Bio.motifs.jaspar", _EXTRA)  # dnachisel's pseudocount routine expects Bio.motifs.jaspar to be loaded
    scan = None if scan_sequence is None else _dna(scan_sequence, "scan_sequence")
    try:
        patterns = dnachisel.MotifPssmPattern.list_from_file(
            str(source), fmt, threshold=t, pseudocounts=_pseudocounts(pseudocounts), relative_threshold=r
        )
    except (ValueError, TypeError, OSError, KeyError) as exc:
        raise ConfigurationError(f"DnaChisel could not read the motif file: {exc}") from exc
    return {
        "status": "completed",
        "operation": "read_pssm_patterns_from_motif_file",
        "file_format": fmt,
        "motif_count": len(patterns),
        "motifs": [_pssm_payload(p, scan) for p in patterns],
        "provenance": h.provenance("dnachisel", "dnachisel.MotifPssmPattern.list_from_file"),
    }


__all__ = [
    "build_pssm_pattern_from_sequences",
    "evaluate_dnachisel_design_rule",
    "read_pssm_patterns_from_motif_file",
]

_register_public_verbose(globals(), tuple(__all__))
