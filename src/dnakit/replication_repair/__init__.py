"""Local candidate adapters for the replication repair capability family."""

from __future__ import annotations

from collections.abc import Mapping

from dnakit._candidate_runtime import run_candidate_operation as run_candidate_operation

from . import _impl_repl as _impl


def annotate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Annotate genomic features as leading- or lagging-strand relative to the replication origin and terminus.

    input: mapping {'features': [{'id', 'start', 'end', 'strand'}], 'length': chromosome length (or 'sequence': DNA), optional
    'origin' and 'terminus'}. Coordinates are 0-based, half-open; the chromosome is treated as circular.
    options: origin, terminus (override the input; when no origin is given it is estimated from the cumulative GC skew of
    'sequence', a heuristic; the terminus defaults to the antipode).
    A feature (midpoint) is leading-strand (co-directional with the fork) when it is on the + strand in the right replichore
    (origin to terminus, increasing coordinates) or on the - strand in the left replichore. Result: per-feature replichore,
    orientation and distance, plus the leading-strand fraction.
    """
    return _impl.run_op("annotate", _impl.op_annotate, input, options)


def build_recombination_map(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Build a genetic (recombination) map: cM distances and cM/Mb rates between ordered markers.

    input: list (or {'markers': [...]}) of at least two markers sorted by position, each {'id', 'position' (bp)} plus, for every
    marker after the first, one of 'cm' (cumulative cM), 'recombination_fraction' (r to the previous marker) or
    'recombinants' with 'meioses' (r = recombinants / meioses).
    options: map_function ('haldane' default: d = -1/2 ln(1-2r); 'kosambi': d = 1/4 ln((1+2r)/(1-2r)); 'none': d = r), applied
    only to r-based markers.
    Result: per-interval delta cM, cumulative cM and rate = delta cM / interval Mb, plus the total and mean rate.
    """
    return _impl.run_op("build_recombination_map", _impl.op_build_recombination_map, input, options)


def calculate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Calculate windowed GC (or AT) skew and cumulative skew of a DNA sequence and a heuristic replication origin/terminus estimate.

    input: DNA string, single-record FASTA text, a file path or {'sequence': ...} (up to 20 Mb).
    options: window (default min(1000, length)), step (default window), skew_type ('GC' = (G-C)/(G+C), default, or 'AT' = (A-T)/(A+T)),
    include_windows (default True), max_windows_returned (default 20000).
    Cumulative skew is the running sum of the window skews. Heuristic: the replication origin is estimated as the centre of the window
    with the minimum cumulative skew and the terminus as the maximum (Lobry 1996, Grigoriev 1998); this only holds for bacterial-like
    chromosomes with a strand composition bias and has window-size resolution.
    """
    return _impl.run_op("calculate", _impl.op_calculate, input, options)


def classify_repair_pathway(
    input: object, *, options: Mapping[str, object] | None = None
) -> object:
    """Classify the likely DSB repair pathway (NHEJ, MMEJ, SSA, HDR) from the reference and repaired sequences with an indel-signature rule.

    input: mapping {'reference': unedited sequence, 'edited': repaired sequence, optional 'donor': expected HDR product}.
    options: min_microhomology (default 2), ssa_min_homology (default 25), max_nhej_insertion (default 20).
    Heuristic rule, not a measurement: edited == donor -> HDR; a pure deletion whose junctions share >= 25 identical bp -> SSA, >= 2 bp
    -> MMEJ, otherwise NHEJ; a pure insertion up to 20 bp -> NHEJ; substitutions, complex indels and long insertions -> 'undetermined'.
    Junction homology is how far the deletion can slide (prefix match + suffix match - edited length).
    """
    return _impl.run_op("classify_repair_pathway", _impl.op_classify_repair_pathway, input, options)


def compare(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Compare the cumulative GC (or AT) skew profiles and origin estimates of two DNA sequences.

    input: {'a': sequence, 'b': sequence} or a list of exactly two DNA strings/FASTA texts.
    options: window (default length/20, at most 1000), skew_type ('GC' default or 'AT'), points (default 100).
    Cumulative skews are resampled to a common number of points over relative position and compared by Pearson correlation
    (-1 for an inverted skew, e.g. a complementary sequence); the relative origin positions (skew minima, a heuristic) are compared.
    """
    return _impl.run_op("compare", _impl.op_compare, input, options)


def detect(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Detect mismatch-repair-relevant motifs in a DNA sequence: GATC Dam sites and perfect mono- to tetranucleotide microsatellite tracts.

    input: DNA string, single-record FASTA text, a file path or {'sequence': ...}.
    options: min_copies (mapping unit length 1-6 -> minimum copies; defaults {1: 8, 2: 5, 3: 4, 4: 3}), max_results (default 1000).
    GATC is the methylation site E. coli MutHLS uses for strand discrimination; microsatellite tracts are the loci that mismatch-repair
    deficient cells expand or contract. This is an exact string scan; it does not model instability.
    """
    return _impl.run_op("detect", _impl.op_detect, input, options)


def map(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Map replication timing along the genome as the median-normalised log2 ratio of S-phase to G1 read depth.

    input: mapping {'s_phase': [depth per window], 'g1': [depth per window]} over the same fixed-size windows.
    options: pseudocount (default 0; required > 0 if windows can be empty), smoothing_windows (odd, default 1 = none; running median),
    window_size (bp, default 1, used for coordinates), center (default 0), max_windows_returned (default 20000).
    Score = log2((S/median S)/(G1/median G1)); higher means earlier replication (Koren et al. 2012 style). Windows above `center` are
    'early', below are 'late'; local maxima above `center` are reported as early-replication peaks.
    """
    return _impl.run_op("map", _impl.op_map, input, options)


def plot(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Plot the cumulative GC (or AT) skew of a DNA sequence with the heuristic origin and terminus marked, as an SVG string.

    input and skew options as calculate (window, step, skew_type). Extra options: output_path (optional .svg file; refuses to
    overwrite unless overwrite=True). The origin/terminus marks are the cumulative-skew extrema (heuristic).
    """
    return _impl.run_op("plot", _impl.op_plot, input, options)


def predict(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Predict candidate bacterial replication-origin (oriC) regions by DnaA-box (TTATNCACA) density.

    input: DNA string, single-record FASTA text, a file path or {'sequence': ...}.
    options: max_mismatches (Hamming vs TTATNCACA, default 0), window (default 500 bp), top_n (default 3).
    Heuristic: the non-overlapping windows with the most DnaA boxes on either strand are returned, together with the GC-skew origin
    estimate and each window's distance to it for cross-checking. Box density alone is a weak, bacteria-specific predictor.
    """
    return _impl.run_op("predict", _impl.op_predict, input, options)


def simulate(input: object, *, options: Mapping[str, object] | None = None) -> object:
    """Simulate replication-fork progression from one or more origins and return the replication-timing profile.

    input: genome length in bp (int, or {'length': n}) or a DNA sequence (its length is used).
    options: origins (list of positions or {'position', 'firing_time'}; default n_origins random positions from seed), n_origins (1),
    seed (0), fork_speed (bp per time unit, default 1000), circular (default False), window (default length/100), time_points (5).
    Deterministic model: replication time of x = min over origins of firing_time + |x - origin| / fork_speed. Result: timing per
    window, S-phase duration, termination sites where forks converge and the fraction replicated over time. No stochastic firing or stalling.
    """
    return _impl.run_op("simulate", _impl.op_simulate, input, options)


__all__ = [
    "annotate",
    "build_recombination_map",
    "calculate",
    "classify_repair_pathway",
    "compare",
    "detect",
    "map",
    "plot",
    "predict",
    "simulate",
]

__candidate_adapter_operations__ = set(
    [
        "annotate",
        "build_recombination_map",
        "calculate",
        "classify_repair_pathway",
        "compare",
        "detect",
        "map",
        "plot",
        "predict",
        "simulate",
    ]
)


from dnakit._verbosity import register_public as _register_public_verbose

_register_public_verbose(globals())
