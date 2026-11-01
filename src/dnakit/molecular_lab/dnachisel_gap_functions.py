"""DnaChisel-era DNA operations that DNAKit did not yet expose (wave 1a coverage gap).

The functions here close the gaps found when the DnaChisel (``dnachisel`` 3.x) rows of the
Fig. 2 coverage worklist were rechecked against DNAKit. Each one either calls the original
DnaChisel or Biopython routine, or implements the same documented definition directly, and
returns a JSON-safe dictionary.

Groups:

* Composition and sequence rules: ``generate_random_dna_with_composition`` (gc_share/probas),
  ``gc_fraction_or_windows`` (DnaChisel's global and sliding-window GC), ``load_dna_record``
  (FASTA, GenBank, EMBL and SnapGene by extension), ``annotate_pattern_occurrences_features``
  (IUPAC pattern occurrences),
  ``list_common_restriction_enzymes_filtered`` (site lengths, optimum temperatures, unlike sites).
* Constraint reports: ``constraints_breaches_table`` (several sequences in one table) and the
  breach plots ``plot_constraint_breaches_png`` and ``write_constraint_breaches_pdf``.
* Optimization reports and plots: ``write_dnachisel_optimization_report``,
  ``write_dnachisel_no_solution_report``, ``plot_optimization_changes_png``. They need the
  optional report packages and return the output path plus the underlying data.
* External programs: ``blast_sequence_hits`` (blastn), ``create_bowtie_index`` and
  ``find_bowtie_matches`` (bowtie-build and bowtie 1). They run through the bounded runner with
  fixed argument lists and never through a shell. A missing program raises
  ``BackendUnavailableError`` with an install hint.

Positions are 0-based; segments and intervals are half-open ``[start, end)``. Sequences are
upper-cased (DnaChisel itself requires upper case for its GC function; see each docstring).
"""

from __future__ import annotations

import csv
import importlib.metadata
import math
import os
import re
import shutil
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.backends.external import execute_bounded_command
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, BackendUnavailableError, ConfigurationError
from dnakit.molecular_lab.backends import _dnachisel_problem, _dnachisel_specs
from dnakit.molecular_lab.dnachisel_functions import find_sequence_difference_segments

_EXTRA = "bio-dnachisel"
_PLOT_EXTRA = "bio-plots"
_DNA_FEATURES_EXTRA = "bio-dna-features-viewer"
_BIO_EXTRA = "bio-formats"
_SNAPGENE_EXTRA = "bio-snapgene"
_GENEBLOCKS_EXTRA = "bio-geneblocks"
_REPORT_EXTRA = "bio-dnachisel-reports"
_MAX_DNA = 200_000
_MAX_SCAN_CELLS = 50_000_000
_MAX_RANDOM_LENGTH = 200_000
_MAX_SEQUENCES = 1_000
_MAX_CONSTRAINT_TOTAL_BASES = 2_000_000
_MAX_FILE_BYTES = 50_000_000
_MAX_BOWTIE_BASES = 10_000_000
_IUPAC = frozenset("ACGTWSMKRYBDHVN")
_IUPAC_CLASS = {
    "A": "A", "C": "C", "G": "G", "T": "T", "W": "AT", "S": "CG", "M": "AC", "K": "GT",
    "R": "AG", "Y": "CT", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
_IUPAC_COMPLEMENT = str.maketrans("ACGTWSMKRYBDHVN", "TGCAWSKMYRVHDBN")
_BLAST_TASKS = ("blastn", "blastn-short", "dc-megablast", "megablast")
_DUST = re.compile(r"^(no|yes|\d+ \d+ \d+)$")
_BREACH = re.compile(r"(\d+)-(\d+)(\(\+\)|\(-\)|)")
_VERSION = re.compile(r"(\d+\.\d+\.\d+)")
_TIMEOUT_SECONDS = 3600.0
_MAX_LOG_CHARS = 4000
_MAX_OUTPUT_BYTES = 50_000_000
_MAX_ARTIFACT_BYTES = 2_000_000_000


def _integer(value: object, name: str, *, minimum: int, maximum: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not minimum <= value <= maximum:
        raise ConfigurationError(f"{name} must be an integer in [{minimum}, {maximum}].")
    return value


def _number(value: object, name: str, *, minimum: float, maximum: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigurationError(f"{name} must be a number in [{minimum}, {maximum}].")
    number = float(value)
    if not math.isfinite(number) or not minimum <= number <= maximum:
        raise ConfigurationError(f"{name} must be a finite number in [{minimum}, {maximum}].")
    return number


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be a boolean.")
    return value


def _iupac_text(value: object, name: str = "sequence", *, maximum: int = _MAX_DNA, allow_empty: bool = False) -> str:
    """Return an upper-case IUPAC DNA string (whitespace removed); reject anything else."""

    if not isinstance(value, str):
        raise ConfigurationError(f"{name} must be a DNA string.")
    text = "".join(value.split()).upper()
    if not text and not allow_empty:
        raise ConfigurationError(f"{name} must be a non-empty DNA string.")
    if len(text) > maximum:
        raise ConfigurationError(f"{name} exceeds {maximum} bases.")
    if set(text) - _IUPAC:
        raise ConfigurationError(f"{name} must contain only IUPAC DNA letters.")
    return text


def _provenance(**notes: object) -> dict[str, Any]:
    try:
        installed = importlib.metadata.version("dnachisel")
    except importlib.metadata.PackageNotFoundError:
        installed = None
    return {"package": "dnachisel", "version": installed, **notes}


def _fresh_path(value: object, name: str) -> Path:
    """Return a path that does not exist yet; its parent is created. Existing outputs are never reused."""

    if not isinstance(value, (str, os.PathLike)) or not str(value).strip():
        raise ConfigurationError(f"{name} must be a path.")
    path = Path(value).expanduser().resolve()
    if path.exists():
        raise ConfigurationError(f"{name} already exists: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def _relative_files(root: Path) -> list[str]:
    return sorted(str(item.relative_to(root)) for item in root.rglob("*") if item.is_file())


def _csv_rows(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            parsed: dict[str, Any] = {}
            for key, value in row.items():
                text = "" if value is None else str(value)
                parsed[key] = int(text) if re.fullmatch(r"-?\d+", text) else text
            rows.append(parsed)
    return rows


def _rc(text: str) -> str:
    return text.translate(_IUPAC_COMPLEMENT)[::-1]


def _named_sequences(value: object) -> list[tuple[str, str]]:
    """Accept {name: sequence} or [{"name": ..., "sequence": ...}] and return unique (name, DNA) pairs."""

    if isinstance(value, Mapping):
        items = list(value.items())
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        items = []
        for entry in value:
            if not isinstance(entry, Mapping) or not isinstance(entry.get("name"), str):
                raise ConfigurationError("Each sequence needs a string 'name' and a 'sequence'.")
            items.append((entry["name"], entry.get("sequence")))
    else:
        raise ConfigurationError("sequences must be a mapping or a list of name/sequence mappings.")
    if not items or len(items) > _MAX_SEQUENCES:
        raise ConfigurationError(f"sequences must contain 1 to {_MAX_SEQUENCES} entries.")
    pairs: list[tuple[str, str]] = []
    names: set[str] = set()
    total = 0
    for name, sequence in items:
        if not isinstance(name, str) or not name.strip() or name in names:
            raise ConfigurationError("Sequence names must be unique non-empty strings.")
        names.add(name)
        text = _iupac_text(sequence, f"sequence {name!r}")
        total += len(text)
        if total > _MAX_CONSTRAINT_TOTAL_BASES:
            raise ConfigurationError(f"sequences exceed {_MAX_CONSTRAINT_TOTAL_BASES} bases in total.")
        pairs.append((name, text))
    return pairs


def _parse_breaches(text: str) -> list[list[object]]:
    locations: list[list[object]] = []
    for match in _BREACH.finditer(text):
        strand = {"(+)": 1, "(-)": -1}.get(match.group(3))
        locations.append([int(match.group(1)), int(match.group(2)), strand])
    return locations


def _external_program(name: str, path: str | None, hint: str) -> str:
    """Resolve an executable on PATH or at an explicit path. Never a shell lookup."""

    if path is None or path == "":
        found = shutil.which(name)
        if found is None:
            raise BackendUnavailableError(
                f"External executable {name!r} was not found on PATH.",
                hint=hint,
                code="DNACHISEL_EXTERNAL_UNAVAILABLE",
            )
        return str(Path(found).absolute())
    candidate = Path(path).expanduser()
    if not candidate.is_file() or not os.access(candidate, os.X_OK):
        raise BackendUnavailableError(
            f"External executable is missing or not executable: {path}",
            hint=hint,
            code="DNACHISEL_EXTERNAL_UNAVAILABLE",
        )
    return str(candidate.absolute())


def _execute(
    program: str,
    executable: str,
    arguments: Sequence[str],
    *,
    cwd: Path,
    outputs: Sequence[Path] = (),
    max_output_bytes: int = _MAX_OUTPUT_BYTES,
    verbose: bool,
) -> str:
    result = execute_bounded_command(
        executable,
        tuple(arguments),
        backend_id=program,
        cwd=cwd,
        timeout_seconds=_TIMEOUT_SECONDS,
        max_output_bytes=max_output_bytes,
        monitored_output_paths=tuple(outputs),
        max_monitored_output_bytes=_MAX_ARTIFACT_BYTES,
        verbose=verbose,
    )
    if result.return_code != 0:
        raise BackendExecutionError(
            f"{program} exited with code {result.return_code}.",
            code="DNACHISEL_EXTERNAL_FAILED",
            context={"log_tail": result.output[-_MAX_LOG_CHARS:]},
        )
    missing = [str(path) for path in outputs if not path.is_file()]
    if missing:
        raise BackendExecutionError(
            f"{program} did not write the expected files: {missing}",
            code="DNACHISEL_EXTERNAL_OUTPUT_MISSING",
        )
    return result.output


def _program_version(program: str, executable: str, arguments: Sequence[str], cwd: Path, verbose: bool) -> str | None:
    output = _execute(program, executable, arguments, cwd=cwd, verbose=verbose)
    match = _VERSION.search(output)
    return match.group(1) if match else None


def _solve(problem: Any, *, has_objectives: bool) -> None:
    """Resolve constraints, then optimise objectives when there are any (as optimize_dnachisel_record)."""

    problem.resolve_constraints()
    if has_objectives:
        problem.optimize()


# ---------------------------------------------------------------------------
# Composition and sequence rules
# ---------------------------------------------------------------------------


def generate_random_dna_with_composition(
    length: int,
    *,
    gc_share: float | None = None,
    probas: Mapping[str, float] | None = None,
    seed: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Generate a random DNA sequence with a GC share or per-base probabilities.

    ``gc_share`` gives the total G+C probability (split evenly between G and C, and A and T).
    ``probas`` maps bases to probabilities (for example ``{"A": 0.2, "T": 0.3, "G": 0.3, "C": 0.2}``)
    and must sum to 1. Without either, all four bases are equiprobable. Giving both is an error.
    DnaChisel draws with numpy's global generator, so ``seed`` seeds that generator (a side effect
    of the original routine, kept for reproducibility).
    """

    dnachisel = require("dnachisel", _EXTRA)
    size = _integer(length, "length", minimum=0, maximum=_MAX_RANDOM_LENGTH)
    share = None if gc_share is None else _number(gc_share, "gc_share", minimum=0.0, maximum=1.0)
    if probas is not None and share is not None:
        raise ConfigurationError("Give gc_share or probas, not both.")
    weights: dict[str, float] | None = None
    if probas is not None:
        if not isinstance(probas, Mapping) or not probas:
            raise ConfigurationError("probas must be a non-empty mapping of base to probability.")
        if set(probas) - set("ACGT") or len(set(probas)) != len(probas):
            raise ConfigurationError("probas keys must be distinct bases among A, C, G and T.")
        weights = {str(key): _number(value, f"probas[{key}]", minimum=0.0, maximum=1.0) for key, value in probas.items()}
        if abs(sum(weights.values()) - 1.0) > 1e-6:
            raise ConfigurationError("probas must sum to 1.")
    seed_value = None if seed is None else _integer(seed, "seed", minimum=0, maximum=2**32 - 1)
    try:
        sequence = str(
            dnachisel.random_dna_sequence(size, gc_share=share, probas=weights, seed=seed_value)
        )
    except (ValueError, TypeError) as exc:
        raise ConfigurationError(f"DnaChisel could not generate the sequence: {exc}") from exc
    gc_count = sequence.count("G") + sequence.count("C")
    source = "gc_share" if share is not None else ("probas" if weights is not None else "uniform")
    return {
        "status": "completed",
        "operation": "generate_random_dna_with_composition",
        "sequence": sequence,
        "length": len(sequence),
        "composition_source": source,
        "gc_share": share,
        "probas": weights,
        "gc_fraction_observed": (gc_count / len(sequence)) if sequence else None,
        "seed": seed_value,
        "provenance": _provenance(method="dnachisel.random_dna_sequence (numpy random.choice)"),
    }


def gc_fraction_or_windows(
    sequence: str,
    *,
    window_size: int | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return the global GC fraction, or the GC fraction of every sliding window.

    Global: (count of G and C) / len(sequence), where every symbol, including N and other IUPAC
    codes, counts in the denominator (DnaChisel's definition). With ``window_size`` the result is
    one value per window start ``i`` (``0 <= i <= len - window_size``), covering ``[i, i + window_size)``.
    The sequence is upper-cased first, because DnaChisel counts only upper-case G and C.
    """

    module = require("dnachisel.biotools.gc_content", _EXTRA)
    text = _iupac_text(sequence, "sequence", allow_empty=window_size is not None)
    if window_size is None:
        fraction = float(module.gc_content(text))
        return {
            "status": "completed",
            "operation": "gc_fraction_or_windows",
            "sequence_length": len(text),
            "window_size": None,
            "gc_fraction": fraction,
            "window_count": 0,
            "window_gc_fractions": [],
            "provenance": _provenance(method="dnachisel.biotools.gc_content.gc_content, global"),
        }
    size = _integer(window_size, "window_size", minimum=1, maximum=_MAX_DNA)
    values: list[float] = []
    if size <= len(text):
        values = [float(item) for item in module.gc_content(text, window_size=size)]
    return {
        "status": "completed",
        "operation": "gc_fraction_or_windows",
        "sequence_length": len(text),
        "window_size": size,
        "gc_fraction": None,
        "window_count": len(values),
        "window_gc_fractions": values,
        "provenance": _provenance(method="dnachisel.biotools.gc_content.gc_content, sliding window"),
    }


def load_dna_record(
    source: str | Path,
    *,
    file_format: str = "auto",
    name: str | None = None,
    linear: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Read one sequence record from FASTA, GenBank, EMBL or SnapGene (.dna) and return it.

    ``file_format="auto"`` chooses by extension: ``.gb``, ``.gbk`` and ``.genbank`` are GenBank,
    ``.fa``, ``.fasta`` and ``.fna`` are FASTA and ``.dna`` is SnapGene (which needs the
    ``snapgene_reader`` package). A file with more than one record is an error (as DnaChisel's
    ``load_record`` with ``SeqIO.read``). The sequence keeps the case used in the file.
    """

    biopython_formats = {"genbank": "genbank", "embl": "embl", "fasta": "fasta"}
    selected = file_format
    if selected not in {"auto", "genbank", "embl", "fasta", "snapgene"}:
        raise ConfigurationError("file_format must be auto, genbank, embl, fasta or snapgene.")
    path = Path(source).expanduser()
    if not path.is_file():
        raise ConfigurationError(f"Input file does not exist: {path}")
    if path.stat().st_size > _MAX_FILE_BYTES:
        raise ConfigurationError(f"Input file exceeds {_MAX_FILE_BYTES} bytes.")
    if selected == "auto":
        suffix = path.suffix.lower()
        if suffix in {".gb", ".gbk", ".genbank"}:
            selected = "genbank"
        elif suffix in {".fa", ".fasta", ".fna"}:
            selected = "fasta"
        elif suffix == ".dna":
            selected = "snapgene"
        else:
            raise ConfigurationError(f"Unknown format for file: {path.name}; pass file_format.")
    linear_flag = _flag(linear, "linear")
    if selected == "snapgene":
        snapgene = require("snapgene_reader", _SNAPGENE_EXTRA)
        record = snapgene.snapgene_file_to_seqrecord(str(path))
    else:
        seqio = require("Bio.SeqIO", _BIO_EXTRA)
        try:
            record = seqio.read(str(path), biopython_formats[selected])
        except ValueError as exc:
            raise ConfigurationError(f"Could not read exactly one {selected} record: {exc}") from exc
    if name is not None:
        if not isinstance(name, str) or not name.strip():
            raise ConfigurationError("name must be a non-empty string.")
        record.id = name
        record.name = name.replace(" ", "_")[:20]
    sequence = str(record.seq)
    if len(sequence) > _MAX_FILE_BYTES:
        raise ConfigurationError("The record is too long for this function.")
    return {
        "status": "completed",
        "operation": "load_dna_record",
        "source": str(path),
        "format": selected,
        "id": str(record.id),
        "name": str(record.name),
        "description": str(record.description),
        "sequence": sequence,
        "length": len(sequence),
        "topology": "linear" if linear_flag else "circular",
        "feature_count": len(record.features),
        "provenance": _provenance(method=f"Bio.SeqIO {selected}" if selected != "snapgene" else "snapgene_reader"),
    }


def annotate_pattern_occurrences_features(
    sequence: str,
    pattern: str,
    *,
    feature_type: str = "misc_feature",
    prefix: str = "!",
    output_path: str | Path | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find every occurrence of an IUPAC DNA pattern on both strands and annotate them.

    Returns one feature per occurrence (``start`` inclusive, ``end`` exclusive, ``strand`` +1 or -1,
    label ``prefix + pattern``). Overlapping occurrences are kept. A palindromic pattern is searched
    on the forward strand only, as DnaChisel does. Optionally a new GenBank file is written at
    ``output_path``. Each IUPAC code is matched by its full set of bases (N is any base, R is A or G,
    and so on). DnaChisel's own pattern regex maps N, R and Y to the wrong sets, so its annotations
    differ for patterns that contain those codes.
    """

    text = _iupac_text(sequence, "sequence", allow_empty=True)
    motif = _iupac_text(pattern, "pattern", maximum=1000)
    if not isinstance(feature_type, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_']{0,39}", feature_type):
        raise ConfigurationError("feature_type must be a GenBank feature key such as misc_feature.")
    if not isinstance(prefix, str) or len(prefix) > 20:
        raise ConfigurationError("prefix must be a string of at most 20 characters.")
    if len(text) * len(motif) > _MAX_SCAN_CELLS:
        raise ConfigurationError("sequence and pattern are too large to scan in one call.")
    label = prefix + motif
    occurrences: list[dict[str, Any]] = []

    def scan(expression: str, strand: int) -> None:
        for match in re.finditer(r"(?=(" + expression + "))", text):
            start = match.start()
            occurrences.append({"type": feature_type, "start": start, "end": start + len(motif), "strand": strand, "label": label})

    forward = "".join(f"[{_IUPAC_CLASS[symbol]}]" for symbol in motif)
    scan(forward, 1)
    reverse_motif = _rc(motif)
    if reverse_motif != motif:
        scan("".join(f"[{_IUPAC_CLASS[symbol]}]" for symbol in reverse_motif), -1)
    occurrences.sort(key=lambda item: (item["start"], item["end"], -item["strand"]))
    result: dict[str, Any] = {
        "status": "completed",
        "operation": "annotate_pattern_occurrences_features",
        "sequence_length": len(text),
        "pattern": motif,
        "palindromic": reverse_motif == motif,
        "occurrence_count": len(occurrences),
        "features": occurrences,
        "output_path": None,
        "provenance": _provenance(method="IUPAC class regex on both strands (overlapping)"),
    }
    if output_path is not None:
        target = _fresh_path(output_path, "output_path")
        if target.suffix.lower() not in {".gb", ".gbk", ".genbank"}:
            raise ConfigurationError("output_path must end with .gb, .gbk or .genbank.")
        seqio = require("Bio.SeqIO", _BIO_EXTRA)
        seq_record = require("Bio.SeqRecord", _BIO_EXTRA)
        seq_feature = require("Bio.SeqFeature", _BIO_EXTRA)
        biopython_seq = require("Bio.Seq", _BIO_EXTRA)
        features = [
            seq_feature.SeqFeature(
                seq_feature.FeatureLocation(item["start"], item["end"], item["strand"]),
                type=feature_type,
                qualifiers={"label": [label]},
            )
            for item in occurrences
        ]
        record = seq_record.SeqRecord(
            biopython_seq.Seq(text),
            id="record",
            name="record",
            description="pattern occurrences",
            features=features,
            annotations={"molecule_type": "DNA"},
        )
        seqio.write(record, str(target), "genbank")
        result["output_path"] = str(target)
    return result


def list_common_restriction_enzymes_filtered(
    *,
    site_length: Sequence[int] = (6,),
    opt_temp: Sequence[int] = (37,),
    min_suppliers: int = 1,
    site_unlike: Sequence[str] = (),
    verbose: bool = True,
) -> dict[str, Any]:
    """List restriction enzymes by site length, optimum temperature, supplier count and unlike sites.

    This is DnaChisel's ``list_common_enzymes`` with its full parameter set (``DNAKit``'s
    ``list_common_restriction_enzymes`` keeps only a single site length and the supplier count).
    ``site_unlike`` names enzymes whose recognition sites (all IUPAC variants) must not occur
    among the returned enzymes' sites. Names are returned sorted.
    """

    dnachisel = require("dnachisel", _EXTRA)
    lengths = _int_items(site_length, "site_length", minimum=4, maximum=12)
    temps = _int_items(opt_temp, "opt_temp", minimum=0, maximum=120)
    suppliers = _integer(min_suppliers, "min_suppliers", minimum=0, maximum=20)
    if isinstance(site_unlike, (str, bytes)) or not isinstance(site_unlike, Sequence):
        raise ConfigurationError("site_unlike must be a list of enzyme names.")
    unlike = []
    for item in site_unlike:
        if not isinstance(item, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]{0,40}", item):
            raise ConfigurationError("Each site_unlike entry must be an enzyme name.")
        unlike.append(item)
    try:
        names = sorted(
            str(name)
            for name in dnachisel.list_common_enzymes(
                site_length=lengths, opt_temp=temps, min_suppliers=suppliers, site_unlike=tuple(unlike)
            )
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise ConfigurationError(f"DnaChisel could not list enzymes: {exc}") from exc
    return {
        "status": "completed",
        "operation": "list_common_restriction_enzymes_filtered",
        "site_length": list(lengths),
        "opt_temp": list(temps),
        "min_suppliers": suppliers,
        "site_unlike": unlike,
        "enzyme_count": len(names),
        "enzymes": names,
        "provenance": _provenance(table="DnaChisel bundled enzyme table (Biopython Restriction)"),
    }


def _int_items(values: object, name: str, *, minimum: int, maximum: int) -> tuple[int, ...]:
    if isinstance(values, (str, bytes)) or not isinstance(values, Sequence) or not values or len(values) > 20:
        raise ConfigurationError(f"{name} must be a list of 1 to 20 integers.")
    return tuple(_integer(item, name, minimum=minimum, maximum=maximum) for item in values)


# ---------------------------------------------------------------------------
# Constraint breach tables and plots
# ---------------------------------------------------------------------------


def _breach_frame(constraints: object, pairs: list[tuple[str, str]], *, display_locations: bool) -> tuple[Any, list[dict[str, Any]], list[str]]:
    dnachisel = require("dnachisel", _EXTRA)
    breaches = require("dnachisel.reports.constraints_reports.constraints_breaches_dataframe", _EXTRA)
    specs = _dnachisel_specs(constraints, "constraints", dnachisel)
    if not specs:
        raise ConfigurationError("Provide at least one constraint.")
    frame = breaches.constraints_breaches_dataframe(specs, pairs, display_constraints_locations=display_locations)
    labels = [str(column) for column in frame.columns]
    rows = []
    for name, row in frame.iterrows():
        entry: dict[str, Any] = {"sequence": str(name)}
        for label in labels:
            entry[label] = "" if row[label] is None else str(row[label])
        rows.append(entry)
    return frame, rows, labels


def constraints_breaches_table(
    sequences: object,
    *,
    constraints: Sequence[Mapping[str, Any]],
    display_constraints_locations: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Return a table of constraint breaches for several sequences (one row per sequence).

    ``sequences`` is a mapping ``{name: sequence}`` or a list of ``{"name", "sequence"}``
    mappings. Each column is one constraint, labelled in DnaChisel's breach form; a cell lists the
    breached locations as ``start-end(strand)`` separated by commas, and is empty when the
    constraint passes. ``constraints`` uses the same allow-list as ``evaluate_dnachisel_specs``.
    """

    pairs = _named_sequences(sequences)
    flag = _flag(display_constraints_locations, "display_constraints_locations")
    _, rows, labels = _breach_frame(constraints, pairs, display_locations=flag)
    return {
        "status": "completed",
        "operation": "constraints_breaches_table",
        "sequence_count": len(rows),
        "constraint_labels": labels,
        "rows": rows,
        "provenance": _provenance(method="dnachisel.reports.constraints_reports.constraints_breaches_dataframe"),
    }


def plot_constraint_breaches_png(
    sequence: str,
    *,
    output_path: str | Path,
    constraints: Sequence[Mapping[str, Any]],
    name: str = "sequence",
    figure_width: float = 10.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw one sequence with its constraint breaches highlighted and save it as a PNG.

    Uses DnaChisel's ``plot_breaches_record`` (matplotlib and dna-features-viewer). The breached
    locations are returned as data next to the image path.
    """

    require("matplotlib", _PLOT_EXTRA)
    plt = require("matplotlib.pyplot", _PLOT_EXTRA)
    require("dna_features_viewer", _DNA_FEATURES_EXTRA)
    breach_plots = require("dnachisel.reports.constraints_reports.constraints_reports", _EXTRA)
    text = _iupac_text(sequence, "sequence")
    width = _number(figure_width, "figure_width", minimum=2.0, maximum=40.0)
    target = _fresh_path(output_path, "output_path")
    if target.suffix.lower() != ".png":
        raise ConfigurationError("output_path must end with .png.")
    pairs = _named_sequences([{"name": name, "sequence": text}])
    frame, rows, labels = _breach_frame(constraints, pairs, display_locations=False)
    records = breach_plots.records_from_breaches_dataframe(frame, pairs)
    axes = breach_plots.plot_breaches_record(records[0], figure_width=width)
    figure = axes.figure
    figure.savefig(str(target), bbox_inches="tight")
    plt.close(figure)
    locations = {label: _parse_breaches(rows[0][label]) for label in labels}
    return {
        "status": "completed",
        "operation": "plot_constraint_breaches_png",
        "image_path": str(target),
        "sequence_length": len(text),
        "breaches": locations,
        "constraint_labels": labels,
        "provenance": _provenance(method="dnachisel.reports.constraints_reports.plot_breaches_record"),
    }


def write_constraint_breaches_pdf(
    sequences: object,
    *,
    output_path: str | Path,
    constraints: Sequence[Mapping[str, Any]],
    figure_width: float = 10.0,
    verbose: bool = True,
) -> dict[str, Any]:
    """Write one breach page per sequence to a new PDF and return the breach table as well.

    Uses DnaChisel's ``breaches_records_to_pdf`` (matplotlib and dna-features-viewer).
    """

    require("matplotlib", _PLOT_EXTRA)
    require("dna_features_viewer", _DNA_FEATURES_EXTRA)
    breach_plots = require("dnachisel.reports.constraints_reports.constraints_reports", _EXTRA)
    pairs = _named_sequences(sequences)
    width = _number(figure_width, "figure_width", minimum=2.0, maximum=40.0)
    target = _fresh_path(output_path, "output_path")
    if target.suffix.lower() != ".pdf":
        raise ConfigurationError("output_path must end with .pdf.")
    frame, rows, labels = _breach_frame(constraints, pairs, display_locations=False)
    records = breach_plots.records_from_breaches_dataframe(frame, pairs)
    breach_plots.breaches_records_to_pdf(records, pdf_path=str(target), figure_width=width, logger=None)
    return {
        "status": "completed",
        "operation": "write_constraint_breaches_pdf",
        "pdf_path": str(target),
        "page_count": len(records),
        "sequence_count": len(pairs),
        "constraint_labels": labels,
        "rows": rows,
        "provenance": _provenance(method="dnachisel.reports.constraints_reports.breaches_records_to_pdf"),
    }


# ---------------------------------------------------------------------------
# Optimization reports and plots (need the optional report packages)
# ---------------------------------------------------------------------------


def write_dnachisel_optimization_report(
    sequence: str,
    *,
    output_dir: str | Path,
    constraints: Sequence[Mapping[str, Any]] = (),
    objectives: Sequence[Mapping[str, Any]] = (),
    circular: bool = False,
    project_name: str = "DnaChisel optimization",
    verbose: bool = True,
) -> dict[str, Any]:
    """Optimise a sequence and write DnaChisel's optimization report folder.

    The folder holds the PDF summary, the before/after constraint and objective tables, the
    GenBank files and logs. It is created new and never overwritten. Needs the DnaChisel report
    packages (``pdf_reports`` and ``sequenticon``, extra ``bio-dnachisel-reports``). If constraints
    cannot be satisfied an error is raised; use ``write_dnachisel_no_solution_report`` for that case.
    """

    dnachisel = require("dnachisel", _EXTRA)
    require("pdf_reports", _REPORT_EXTRA)
    require("sequenticon", _REPORT_EXTRA)
    reports = require("dnachisel.reports.optimization_reports", _EXTRA)
    text = _iupac_text(sequence, "sequence")
    specs_c = _dnachisel_specs(constraints, "constraints", dnachisel)
    specs_o = _dnachisel_specs(objectives, "objectives", dnachisel)
    if not specs_c and not specs_o:
        raise ConfigurationError("Provide at least one constraint or objective.")
    flag = _flag(circular, "circular")
    if not isinstance(project_name, str) or not project_name.strip() or len(project_name) > 120:
        raise ConfigurationError("project_name must be a non-empty string of at most 120 characters.")
    target = _fresh_path(output_dir, "output_dir")
    problem = _dnachisel_problem(dnachisel, text, specs_c, specs_o, flag)
    try:
        _solve(problem, has_objectives=bool(specs_o))
    except dnachisel.NoSolutionError as exc:
        raise ConfigurationError(f"No solution: {exc}. Use write_dnachisel_no_solution_report.") from exc
    reports.write_optimization_report(str(target), problem, project_name=project_name)
    return {
        "status": "completed",
        "operation": "write_dnachisel_optimization_report",
        "report_path": str(target),
        "report_files": _relative_files(target),
        "sequence_length": len(problem.sequence),
        "edits": int(problem.number_of_edits()),
        "all_constraints_pass": bool(problem.all_constraints_pass()),
        "constraints_before_after": _csv_rows(target / "constraints_before_and_after.csv"),
        "objectives_before_after": _csv_rows(target / "objectives_before_and_after.csv"),
        "circular": flag,
        "provenance": _provenance(method="dnachisel.reports.optimization_reports.write_optimization_report"),
    }


def write_dnachisel_no_solution_report(
    sequence: str,
    *,
    output_dir: str | Path,
    constraints: Sequence[Mapping[str, Any]] = (),
    objectives: Sequence[Mapping[str, Any]] = (),
    circular: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Report where constraints cannot be satisfied, or report that they can.

    If DnaChisel finds no solution it writes the no-solution report folder (``plots.pdf``, the local
    breach GenBank file and logs) to a new ``output_dir`` and returns the breach location. If the
    constraints are satisfiable the status is ``solved`` and nothing is written.
    """

    dnachisel = require("dnachisel", _EXTRA)
    require("matplotlib", _PLOT_EXTRA)
    require("dna_features_viewer", _DNA_FEATURES_EXTRA)
    reports = require("dnachisel.reports.optimization_reports", _EXTRA)
    text = _iupac_text(sequence, "sequence")
    specs_c = _dnachisel_specs(constraints, "constraints", dnachisel)
    specs_o = _dnachisel_specs(objectives, "objectives", dnachisel)
    if not specs_c and not specs_o:
        raise ConfigurationError("Provide at least one constraint or objective.")
    flag = _flag(circular, "circular")
    target = _fresh_path(output_dir, "output_dir")
    problem = _dnachisel_problem(dnachisel, text, specs_c, specs_o, flag)
    try:
        _solve(problem, has_objectives=bool(specs_o))
    except dnachisel.NoSolutionError as error:
        location = getattr(error, "location", None)
        if location is None:
            return {
                "status": "no_solution_without_location",
                "operation": "write_dnachisel_no_solution_report",
                "report_path": None,
                "report_files": [],
                "message": str(error),
                "location": None,
                "provenance": _provenance(method="dnachisel.NoSolutionError"),
            }
        reports.write_no_solution_report(str(target), problem, error)
        start, end, strand = location.to_tuple()
        return {
            "status": "no_solution",
            "operation": "write_dnachisel_no_solution_report",
            "report_path": str(target),
            "report_files": _relative_files(target),
            "message": str(error),
            "location": [int(start), int(end), strand],
            "circular": flag,
            "provenance": _provenance(method="dnachisel.reports.optimization_reports.write_no_solution_report"),
        }
    return {
        "status": "solved",
        "operation": "write_dnachisel_no_solution_report",
        "report_path": None,
        "report_files": [],
        "message": "The constraints can be satisfied; no report was written.",
        "location": None,
        "circular": flag,
        "provenance": _provenance(method="dnachisel.DnaOptimizationProblem.resolve_constraints"),
    }


def plot_optimization_changes_png(
    sequence: str,
    *,
    output_path: str | Path,
    constraints: Sequence[Mapping[str, Any]] = (),
    objectives: Sequence[Mapping[str, Any]] = (),
    circular: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Optimise a sequence and save a PNG of the changed regions (DnaChisel's plot_optimization_changes).

    Returns the image path, the changed segments (``[start, end)``) and the constraint status.
    Needs the optional ``geneblocks`` package (extra ``bio-geneblocks``).
    """

    dnachisel = require("dnachisel", _EXTRA)
    require("geneblocks", _GENEBLOCKS_EXTRA)
    require("dna_features_viewer", _DNA_FEATURES_EXTRA)
    plt = require("matplotlib.pyplot", _PLOT_EXTRA)
    reports = require("dnachisel.reports.optimization_reports", _EXTRA)
    text = _iupac_text(sequence, "sequence")
    specs_c = _dnachisel_specs(constraints, "constraints", dnachisel)
    specs_o = _dnachisel_specs(objectives, "objectives", dnachisel)
    if not specs_c and not specs_o:
        raise ConfigurationError("Provide at least one constraint or objective.")
    flag = _flag(circular, "circular")
    target = _fresh_path(output_path, "output_path")
    if target.suffix.lower() != ".png":
        raise ConfigurationError("output_path must end with .png.")
    problem = _dnachisel_problem(dnachisel, text, specs_c, specs_o, flag)
    try:
        _solve(problem, has_objectives=bool(specs_o))
    except dnachisel.NoSolutionError as exc:
        raise ConfigurationError(f"No solution: {exc}.") from exc
    axes = reports.plot_optimization_changes(problem)
    figure = axes.figure
    figure.savefig(str(target), bbox_inches="tight")
    plt.close(figure)
    segments = find_sequence_difference_segments(str(problem.sequence_before), str(problem.sequence))
    return {
        "status": "completed",
        "operation": "plot_optimization_changes_png",
        "image_path": str(target),
        "sequence_length": len(problem.sequence),
        "changed_segments": segments["segments"],
        "edits": int(problem.number_of_edits()),
        "all_constraints_pass": bool(problem.all_constraints_pass()),
        "provenance": _provenance(method="dnachisel.reports.optimization_reports.plot_optimization_changes"),
    }


# ---------------------------------------------------------------------------
# External programs: BLAST+ and Bowtie 1
# ---------------------------------------------------------------------------

_BLAST_HINT = (
    "Install the BLAST+ command-line tools (for example `conda install -c bioconda blast`, or the NCBI "
    "release from https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html), then put blastn on "
    "PATH or pass blastn_path. A BLAST database is needed for blast_db (makeblastdb)."
)
_BOWTIE_HINT = (
    "Install Bowtie 1 (for example `conda install -c bioconda bowtie`, or the release from "
    "https://github.com/BenLangmead/bowtie/releases), then put bowtie and bowtie-build on PATH or pass "
    "bowtie_path / bowtie_build_path."
)


def blast_sequence_hits(
    sequence: str,
    *,
    blast_db: str | None = None,
    subject_fasta: str | Path | None = None,
    word_size: int = 4,
    perc_identity: float = 80.0,
    num_alignments: int = 1000,
    ungapped: bool = False,
    num_threads: int = 3,
    culling_limit: int | None = None,
    e_value: float | None = None,
    task: str | None = None,
    dust: str = "no",
    blastn_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """BLAST one sequence against a BLAST database (or a subject FASTA file) with blastn.

    The defaults are DnaChisel's (word size 4, 80 % identity, 1000 alignments, dust off). Output is
    BLAST's XML, parsed into alignments and high-scoring pairs. Exactly one of ``blast_db`` (a database
    prefix) or ``subject_fasta`` must be given. Requires the BLAST+ ``blastn`` program.
    """

    seq_io = require("Bio.Blast.NCBIXML", _BIO_EXTRA)
    text = _iupac_text(sequence, "sequence")
    if (blast_db is None) == (subject_fasta is None):
        raise ConfigurationError("Give exactly one of blast_db or subject_fasta.")
    word = _integer(word_size, "word_size", minimum=4, maximum=28)
    identity = _number(perc_identity, "perc_identity", minimum=0.0, maximum=100.0)
    hits = _integer(num_alignments, "num_alignments", minimum=1, maximum=10_000)
    threads = _integer(num_threads, "num_threads", minimum=1, maximum=64)
    gapless = _flag(ungapped, "ungapped")
    culling = None if culling_limit is None else _integer(culling_limit, "culling_limit", minimum=1, maximum=10_000)
    evalue = None if e_value is None else _number(e_value, "e_value", minimum=0.0, maximum=1e9)
    if task is not None and task not in _BLAST_TASKS:
        raise ConfigurationError(f"task must be one of {', '.join(_BLAST_TASKS)}.")
    if not isinstance(dust, str) or not _DUST.match(dust):
        raise ConfigurationError("dust must be 'no', 'yes' or three integers such as '20 64 1'.")
    database: str | None = None
    subject: Path | None = None
    if blast_db is not None:
        if not isinstance(blast_db, str) or not blast_db.strip() or blast_db.startswith("-"):
            raise ConfigurationError("blast_db must be a database path that does not start with '-'.")
        database = blast_db
    else:
        subject = Path(subject_fasta).expanduser().resolve() if subject_fasta is not None else None
        if subject is None or not subject.is_file():
            raise ConfigurationError("subject_fasta must be an existing FASTA file.")
    executable = _external_program("blastn", blastn_path, _BLAST_HINT)
    with tempfile.TemporaryDirectory(prefix="dnakit_blast_") as directory:
        work = Path(directory)
        query = work / "query.fa"
        query.write_text(f">query\n{text}\n", encoding="utf-8")
        out = work / "blast.xml"
        arguments = [
            "-query", str(query),
            "-out", str(out),
            "-outfmt", "5",
            "-max_target_seqs", str(hits),
            "-word_size", str(word),
            "-num_threads", str(threads),
            "-perc_identity", str(identity),
            "-dust", dust,
        ]
        if database is not None:
            arguments += ["-db", database]
        else:
            arguments += ["-subject", str(subject)]
        if gapless:
            arguments.append("-ungapped")
        if evalue is not None:
            arguments += ["-evalue", str(evalue)]
        if culling is not None:
            arguments += ["-culling_limit", str(culling)]
        if task is not None:
            arguments += ["-task", task]
        version = _program_version("blastn", executable, ["-version"], work, verbose)
        _execute("blastn", executable, arguments, cwd=work, outputs=[out], verbose=verbose)
        with out.open(encoding="utf-8") as handle:
            records = list(seq_io.parse(handle))
        parsed = []
        for record in records:
            alignments = []
            for alignment in record.alignments:
                alignments.append({
                    "title": str(alignment.title),
                    "length": int(alignment.length),
                    "hsps": [
                        {
                            "identities": int(item.identities),
                            "align_length": int(item.align_length),
                            "gaps": int(item.gaps or 0),
                            "query_start": int(item.query_start),
                            "query_end": int(item.query_end),
                            "subject_start": int(item.sbjct_start),
                            "subject_end": int(item.sbjct_end),
                            "score": float(item.score),
                            "evalue": float(item.expect),
                        }
                        for item in alignment.hsps
                    ],
                })
            parsed.append({"query_id": str(record.query), "query_length": int(record.query_length), "alignments": alignments})
    return {
        "status": "completed",
        "operation": "blast_sequence_hits",
        "backend": "blastn",
        "version": version,
        "blast_db": database,
        "subject_fasta": None if subject is None else str(subject),
        "queries": parsed,
        "hit_count": sum(len(item["alignments"]) for item in parsed),
        "provenance": _provenance(method="blastn -outfmt 5 parsed by Bio.Blast.NCBIXML", executable=executable),
    }


def create_bowtie_index(
    sequences: Sequence[str],
    *,
    index_dir: str | Path,
    bowtie_build_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Build a Bowtie 1 index from sequences; sequences are named by their position (0, 1, ...).

    The index is written into a new directory ``index_dir`` as ``bowtie.*.ebwt`` files; the returned
    ``index_prefix`` is what ``find_bowtie_matches`` takes. Requires ``bowtie-build``.
    """

    if isinstance(sequences, (str, bytes)) or not isinstance(sequences, Sequence) or not sequences:
        raise ConfigurationError("sequences must be a non-empty list of DNA strings.")
    if len(sequences) > _MAX_SEQUENCES:
        raise ConfigurationError(f"sequences accepts at most {_MAX_SEQUENCES} entries.")
    texts = [_iupac_text(item, f"sequences[{index}]") for index, item in enumerate(sequences)]
    if sum(len(item) for item in texts) > _MAX_BOWTIE_BASES:
        raise ConfigurationError(f"sequences exceed {_MAX_BOWTIE_BASES} bases in total.")
    target = _fresh_path(index_dir, "index_dir")
    executable = _external_program("bowtie-build", bowtie_build_path, _BOWTIE_HINT)
    target.mkdir(parents=True, exist_ok=False)
    fasta = target / "sequences.fa"
    fasta.write_text("".join(f">{index}\n{item}\n" for index, item in enumerate(texts)), encoding="utf-8")
    prefix = target / "bowtie"
    version = _program_version("bowtie-build", executable, ["--version"], target, verbose)
    _execute(
        "bowtie-build",
        executable,
        ["-f", str(fasta), str(prefix), "--quiet"],
        cwd=target,
        outputs=[target / "bowtie.1.ebwt"],
        verbose=verbose,
    )
    return {
        "status": "completed",
        "operation": "create_bowtie_index",
        "backend": "bowtie-build",
        "version": version,
        "index_prefix": str(prefix),
        "index_files": _relative_files(target),
        "sequence_count": len(texts),
        "provenance": _provenance(method="bowtie-build -f", executable=executable),
    }


def find_bowtie_matches(
    sequence: str,
    *,
    index_prefix: str | Path,
    match_length: int,
    max_mismatches: int = 0,
    bowtie_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Find, for every k-mer of ``sequence``, the best Bowtie 1 hits in an index (DnaChisel's matcher).

    Each k-mer (length ``match_length``) is aligned with at most ``max_mismatches`` mismatches (0 to 3),
    reporting the best alignment only (``--best -k 1``). A match is returned as
    ``{"start", "end", "mismatches"}`` where ``start``/``end`` is the k-mer position in ``sequence``.
    Requires the Bowtie 1 ``bowtie`` program and an index built by ``create_bowtie_index``.
    """

    text = _iupac_text(sequence, "sequence")
    k = _integer(match_length, "match_length", minimum=1, maximum=max(1, len(text)))
    mismatches = _integer(max_mismatches, "max_mismatches", minimum=0, maximum=3)
    prefix = Path(index_prefix).expanduser()
    if not Path(f"{prefix}.1.ebwt").is_file():
        raise ConfigurationError(f"No Bowtie index found at prefix {prefix}; run create_bowtie_index first.")
    kmers = [text[i:i + k] for i in range(len(text) - k + 1)]
    if not kmers:
        return {
            "status": "completed",
            "operation": "find_bowtie_matches",
            "backend": "bowtie",
            "version": None,
            "match_count": 0,
            "matches": [],
            "provenance": _provenance(method="bowtie --best -k 1"),
        }
    executable = _external_program("bowtie", bowtie_path, _BOWTIE_HINT)
    with tempfile.TemporaryDirectory(prefix="dnakit_bowtie_") as directory:
        work = Path(directory)
        reads = work / "kmers.fa"
        reads.write_text("".join(f">{index}\n{kmer}\n" for index, kmer in enumerate(kmers)), encoding="utf-8")
        version = _program_version("bowtie", executable, ["--version"], work, verbose)
        output = _execute(
            "bowtie",
            executable,
            [
                "--best", "-k", "1", "-v", str(mismatches),
                str(prefix), "--quiet", "--suppress", "2,3,4,5,6,7",
                "-f", str(reads),
            ],
            cwd=work,
            verbose=verbose,
        )
    matches: list[dict[str, int]] = []
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) != 2:
            continue
        index = int(fields[0])
        matches.append({"start": index, "end": index + k, "mismatches": fields[1].count(":")})
    matches.sort(key=lambda item: (item["start"], item["mismatches"]))
    return {
        "status": "completed",
        "operation": "find_bowtie_matches",
        "backend": "bowtie",
        "version": version,
        "kmer_length": k,
        "match_count": len(matches),
        "matches": matches,
        "provenance": _provenance(method="bowtie --best -k 1 -v <mismatches> on k-mers", executable=executable),
    }


__all__ = [
    "annotate_pattern_occurrences_features",
    "blast_sequence_hits",
    "constraints_breaches_table",
    "create_bowtie_index",
    "find_bowtie_matches",
    "gc_fraction_or_windows",
    "generate_random_dna_with_composition",
    "list_common_restriction_enzymes_filtered",
    "load_dna_record",
    "plot_constraint_breaches_png",
    "plot_optimization_changes_png",
    "write_constraint_breaches_pdf",
    "write_dnachisel_no_solution_report",
    "write_dnachisel_optimization_report",
]

_register_public_verbose(globals(), tuple(__all__))
