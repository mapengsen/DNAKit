"""Fixed ViennaRNA operations through Biotite 2's typed application interfaces.

No executable name or method is accepted from callers. Optional dependencies are
loaded on use. The upstream default temperature is 37 C; local-fold defaults are
150 bases (RNALfold) and 70 columns (RNALalifold). Results include actual structures,
energies and zero-based spans. The external ViennaRNA programs must be installed.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from importlib.metadata import PackageNotFoundError, version
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as cli
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError, DNAKitError

_EXTRA = "bio-biotite"
_HINT = "Install ViennaRNA (conda install -c bioconda viennarna) and dnakit[bio-biotite]."
_APPS = {
    "fold": ("RNAfold", "FoldApp"),
    "pkplex": ("RNAPKplex", "PKplexApp"),
    "local": ("RNALfold", "LfoldApp"),
    "cofold": ("RNAcofold", "CofoldApp"),
    "duplex": ("RNAduplex", "DuplexApp"),
    "plex": ("RNAplex", "PlexApp"),
    "multifold": ("RNAmultifold", "MultifoldApp"),
    "alifold": ("RNAalifold", "AlifoldApp"),
    "local_alifold": ("RNALalifold", "LalifoldApp"),
    "aliduplex": ("RNAaliduplex", "AliduplexApp"),
}


def _number(value: object, name: str, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a finite number.")
    if not math.isfinite(value) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be between {minimum} and {maximum}.")
    return float(value)


def _strings(values: Sequence[str], *, alignment: bool = False) -> list[str]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not 1 <= len(values) <= 100:
        raise ConfigurationError("Provide between 1 and 100 nucleotide strings.")
    result = []
    for value in values:
        if not isinstance(value, str) or not value or len(value) > 20_000:
            raise ConfigurationError("Each nucleotide string must contain 1 to 20000 bases.")
        text = value.upper().replace("U", "T")
        if set(text) - set("ACGTN-" if alignment else "ACGTN"):
            raise ConfigurationError("Sequences must use A/C/G/T/U/N; only alignments may contain '-'.")
        result.append(text)
    if alignment and (len({len(s) for s in result}) != 1 or any(set(s) == {"-"} for s in result)):
        raise ConfigurationError("Aligned rows must have equal length and contain at least one base.")
    return result


def _plain(result: Any) -> dict[str, Any]:
    record: dict[str, Any] = {"dot_bracket": result.dot_bracket}
    for name in ("free_energy", "covariance_energy"):
        value = getattr(result, name, None)
        record[name] = None if value is None else float(value)
        if value is not None and not math.isfinite(record[name]):
            raise ConfigurationError(f"ViennaRNA returned a non-finite {name}.")
    for name in ("span", "spans"):
        value = getattr(result, name, None)
        if value is not None:
            record[name] = [list(span) for span in value] if name == "spans" else list(value)
    record["base_pairs"] = result.base_pairs().tolist()
    return record


def _run(
    mode: str, groups: Sequence[Sequence[str]], *, temperature: float,
    executable_path: str | None, timeout_seconds: float, max_span: int | None = None,
) -> dict[str, Any]:
    aligned = mode in ("alifold", "local_alifold", "aliduplex")
    texts = [_strings(group, alignment=aligned) for group in groups]
    if aligned and any(len(group) < 2 for group in texts):
        raise ConfigurationError("Each alignment must contain at least two nucleotide rows.")
    temp = _number(temperature, "temperature", -50, 100)
    timeout = _number(timeout_seconds, "timeout_seconds", 0.01, 86400)
    if max_span is not None:
        cli.integer(max_span, "max_span", minimum=1, maximum=20_000)
    program, class_name = _APPS[mode]
    executable = cli.program(program, explicit=executable_path, hint=_HINT)
    sequence = require("biotite.sequence", _EXTRA)
    apps = require("biotite.application.viennarna", _EXTRA)
    if not hasattr(apps, class_name) or not hasattr(getattr(apps, class_name), "run"):
        raise ConfigurationError("These wrappers require Biotite >=2,<3.")
    if aligned:
        align = require("biotite.sequence.align", _EXTRA)
        inputs = [align.Alignment.from_strings(group, sequence.NucleotideSequence) for group in texts]
    elif mode == "multifold":
        inputs = [[sequence.NucleotideSequence(s) for s in texts[0]]]
    else:
        inputs = [sequence.NucleotideSequence(group[0]) for group in texts]
    options: dict[str, Any] = {"temperature": temp}
    if max_span is not None:
        options["max_span"] = max_span
    future = getattr(apps, class_name)(path=executable).run(*inputs, **options)
    try:
        result = future.result(timeout=timeout)
    except TimeoutError:
        future.cancel()
        raise
    except (IndexError, UnicodeError, ValueError) as error:
        raise DNAKitError(
            f"{program} returned output that Biotite could not parse: {error}"
        ) from error
    structures = result if isinstance(result, list) else [result]
    try:
        installed = version("biotite")
    except PackageNotFoundError:
        installed = None
    return {
        "status": "completed", "backend": program, "temperature": temp,
        "structures": [_plain(item) for item in structures],
        "inputs": texts, "span_convention": "zero-based, end-exclusive",
        "provenance": {"package": "biotite", "package_version": installed,
                       "executable": executable, "version": cli._version_of(executable)},
    }


def fold_rna(sequence: str, *, temperature: float = 37, executable_path: str | None = None,
             timeout_seconds: float = 300, verbose: bool = True) -> dict[str, Any]:
    """RNAfold minimum-free-energy structure of one nucleotide strand."""
    return _run("fold", [[sequence]], temperature=temperature, executable_path=executable_path, timeout_seconds=timeout_seconds)


def fold_rna_pseudoknots(sequence: str, *, temperature: float = 37, executable_path: str | None = None,
                        timeout_seconds: float = 300, verbose: bool = True) -> dict[str, Any]:
    """RNAPKplex secondary structure, including pseudoknots."""
    return _run("pkplex", [[sequence]], temperature=temperature, executable_path=executable_path, timeout_seconds=timeout_seconds)


def fold_rna_local(sequence: str, *, max_span: int = 150, temperature: float = 37,
                   executable_path: str | None = None, timeout_seconds: float = 300,
                   verbose: bool = True) -> dict[str, Any]:
    """RNALfold locally stable structures with positions in the original sequence."""
    return _run("local", [[sequence]], max_span=max_span, temperature=temperature, executable_path=executable_path, timeout_seconds=timeout_seconds)


def cofold_rna(sequence1: str, sequence2: str, *, temperature: float = 37,
               executable_path: str | None = None, timeout_seconds: float = 300,
               verbose: bool = True) -> dict[str, Any]:
    """RNAcofold minimum-free-energy structure of two strands, including intramolecular pairs."""
    return _run("cofold", [[sequence1], [sequence2]], temperature=temperature, executable_path=executable_path, timeout_seconds=timeout_seconds)


def duplex_rna(sequence1: str, sequence2: str, *, temperature: float = 37,
               executable_path: str | None = None, timeout_seconds: float = 300,
               verbose: bool = True) -> dict[str, Any]:
    """RNAduplex hybridization restricted to intermolecular base pairs."""
    return _run("duplex", [[sequence1], [sequence2]], temperature=temperature, executable_path=executable_path, timeout_seconds=timeout_seconds)


def find_rna_interactions(target: str, query: str, *, temperature: float = 37,
                         executable_path: str | None = None, timeout_seconds: float = 300,
                         verbose: bool = True) -> dict[str, Any]:
    """RNAplex interaction-site search, returning all reported interactions and spans."""
    return _run("plex", [[target], [query]], temperature=temperature, executable_path=executable_path, timeout_seconds=timeout_seconds)


def fold_rna_complex(sequences: Sequence[str], *, temperature: float = 37,
                     executable_path: str | None = None, timeout_seconds: float = 300,
                     verbose: bool = True) -> dict[str, Any]:
    """RNAmultifold minimum-free-energy structure of a multi-strand complex."""
    if len(_strings(sequences)) < 2:
        raise ConfigurationError("A complex requires at least two strands.")
    return _run("multifold", [sequences], temperature=temperature, executable_path=executable_path, timeout_seconds=timeout_seconds)


def fold_rna_alignment(alignment: Sequence[str], *, temperature: float = 37,
                       executable_path: str | None = None, timeout_seconds: float = 300,
                       verbose: bool = True) -> dict[str, Any]:
    """RNAalifold consensus structure with physical and covariance energies."""
    return _run("alifold", [alignment], temperature=temperature, executable_path=executable_path, timeout_seconds=timeout_seconds)


def fold_rna_alignment_local(alignment: Sequence[str], *, max_span: int = 70,
                             temperature: float = 37, executable_path: str | None = None,
                             timeout_seconds: float = 300, verbose: bool = True) -> dict[str, Any]:
    """RNALalifold locally stable consensus structures and alignment-column spans."""
    return _run("local_alifold", [alignment], max_span=max_span, temperature=temperature, executable_path=executable_path, timeout_seconds=timeout_seconds)


def duplex_rna_alignments(alignment1: Sequence[str], alignment2: Sequence[str], *,
                          temperature: float = 37, executable_path: str | None = None,
                          timeout_seconds: float = 300, verbose: bool = True) -> dict[str, Any]:
    """RNAaliduplex hybridization structure of two nucleotide alignments."""
    return _run("aliduplex", [alignment1, alignment2], temperature=temperature, executable_path=executable_path, timeout_seconds=timeout_seconds)


def rna_structure_coordinates(dot_bracket: str, *, layout: str = "naview",
                               executable_path: str | None = None, timeout_seconds: float = 300,
                               verbose: bool = True) -> dict[str, Any]:
    """RNAplot 2D coordinates for each base of a dot-bracket secondary structure."""
    if not isinstance(dot_bracket, str) or not 1 <= len(dot_bracket) <= 20_000:
        raise ConfigurationError("dot_bracket must contain 1 to 20000 characters.")
    if set(dot_bracket) - set(".()[]{}<>"):
        raise ConfigurationError("dot_bracket contains invalid characters.")
    way = cli.choice(layout, ("radial", "naview", "circular", "rnaturtle", "rnapuzzler"), "layout")
    timeout = _number(timeout_seconds, "timeout_seconds", 0.01, 86400)
    executable = cli.program("RNAplot", explicit=executable_path, hint=_HINT)
    structure = require("biotite.structure", _EXTRA)
    apps = require("biotite.application.viennarna", _EXTRA)
    # Validate bracket pairing before starting an external process.
    structure.base_pairs_from_dot_bracket(dot_bracket)
    future = apps.PlotApp(path=executable).run(dot_bracket, layout_type=apps.PlotApp.Layout[way.upper()])
    try:
        result = future.result(timeout=timeout)
    except TimeoutError:
        future.cancel()
        raise
    # NumPy loadtxt used by Biotite collapses a single coordinate row.
    if getattr(result, "ndim", 2) == 1 and len(dot_bracket) == 1:
        result = result.reshape(1, 2)
    positions = result.tolist()
    if (len(positions) != len(dot_bracket) or
            any(not isinstance(row, list) or len(row) != 2 for row in positions)):
        raise ConfigurationError("RNAplot returned an unexpected number of base coordinates.")
    if any(not math.isfinite(value) for row in positions for value in row):
        raise DNAKitError("RNAplot returned non-finite coordinates for this structure/layout.")
    return {"status": "completed", "backend": "RNAplot", "layout": way,
            "dot_bracket": dot_bracket, "coordinates": positions,
            "provenance": {"executable": executable, "version": cli._version_of(executable)}}


__all__ = ["fold_rna", "fold_rna_pseudoknots", "fold_rna_local", "cofold_rna", "duplex_rna",
           "find_rna_interactions", "fold_rna_complex", "fold_rna_alignment",
           "fold_rna_alignment_local", "duplex_rna_alignments", "rna_structure_coordinates"]
_register_public_verbose(globals(), tuple(__all__))
