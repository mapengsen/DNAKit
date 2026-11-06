"""FIMO (MEME suite) motif scanning with position-specific priors as a typed DNAKit facade.

The user installs the MEME suite (``conda install -c bioconda meme``, or a release from https://meme-suite.org); the
function runs ``fimo`` with an explicit typed signature mapped to fixed flags. No shell is used, the program is resolved
from an explicit path, then PATH, and a missing program raises ``BackendUnavailableError`` with an install hint. The
result is a JSON-serialisable dictionary with status, backend, version, command, return code, output files, the parsed
``fimo.tsv`` rows and a bounded log tail.

With a position-specific prior (``--psp``) and its binned distribution (``--prior-dist``) FIMO scores each occurrence with
the log-posterior odds of Cuellar-Partida et al. 2012 (Bioinformatics 28:56-62) instead of the log-odds score, and
reports p-values and q-values for that score. Defaults are those of the FIMO command line (p-value threshold 1e-4, both
strands, motif pseudocount 0.1, prior ``alpha`` 1.0, 100,000 stored scores). The priors may be a wiggle file (positions
without a value get the median prior) or a MEME PSP file; the distribution file is produced by ``compute-prior-dist`` or
``create-priors`` of the MEME suite.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _gap_cli as _cli
from dnakit.exceptions import ConfigurationError

_HINT = (
    "Install the MEME suite 5.5 or later (conda install -c bioconda meme, or a release from "
    "https://meme-suite.org) and put fimo on PATH or pass fimo_path."
)
_BACKGROUND_SOURCES = {"motif": "--motif--", "uniform": "--uniform--", "nrdb": "--nrdb--"}
_OUTPUT_FILES = (
    "fimo.tsv",
    "fimo.gff",
    "fimo.html",
    "fimo.xml",
    "cisml.xml",
    "best_site.narrowPeak",
)
_TSV_COLUMNS = (
    "motif_id",
    "motif_alt_id",
    "sequence_name",
    "start",
    "stop",
    "strand",
    "score",
    "p_value",
    "q_value",
    "matched_sequence",
)
_INTEGER_COLUMNS = ("start", "stop")
_FLOAT_COLUMNS = ("score", "p_value", "q_value")


def _float_or_none(value: str) -> float | None:
    try:
        return float(value)
    except ValueError:
        return None


def _parse_tsv(path: Path, limit: int) -> tuple[int, list[dict[str, Any]]]:
    """Rows of ``fimo.tsv`` (header and ``#`` comment lines skipped); returns the total count and the first ``limit``."""
    total = 0
    rows: list[dict[str, Any]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("motif_id"):
                continue
            fields = line.split("\t")
            total += 1
            if len(rows) >= limit:
                continue
            row: dict[str, Any] = dict(zip(_TSV_COLUMNS, fields, strict=False))
            for column in _INTEGER_COLUMNS:
                if column in row and row[column].lstrip("-").isdigit():
                    row[column] = int(row[column])
            for column in _FLOAT_COLUMNS:
                if column in row:
                    row[column] = _float_or_none(row[column])
            rows.append(row)
    return total, rows


def _motif_ids(values: Sequence[str] | None) -> list[str]:
    if values is None:
        return []
    if isinstance(values, str) or not values:
        raise ConfigurationError("motif_ids must be a non-empty list of motif identifiers or None.")
    arguments: list[str] = []
    for value in values:
        arguments += ["--motif", _cli.text(value, "motif_ids")]
    return arguments


def scan_motifs_with_position_priors(
    motif_file: str,
    sequence_file: str,
    prior_file: str,
    prior_distribution_file: str,
    output_dir: str,
    *,
    prior_alpha: float = 1.0,
    p_value_threshold: float = 1e-4,
    use_q_value_threshold: bool = False,
    compute_q_values: bool = True,
    scan_reverse_complement: bool = True,
    max_strand: bool = False,
    motif_ids: Sequence[str] | None = None,
    motif_pseudocount: float = 0.1,
    background_file: str | None = None,
    background_source: str | None = None,
    parse_genomic_coordinates: bool = True,
    max_stored_scores: int = 100000,
    max_rows: int = 1000,
    fimo_path: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan DNA sequences for motif occurrences with FIMO using position-specific priors (log-posterior odds).

    Runs ``fimo --psp <prior_file> --prior-dist <prior_distribution_file> [--alpha ...] <motif_file> <sequence_file>``.
    ``motif_file`` holds MEME-format motifs (MEME, STREME, DREME or minimal MEME format), ``sequence_file`` is FASTA.
    ``prior_file`` is a position-specific prior in wiggle or MEME PSP format and ``prior_distribution_file`` is its
    binned distribution (``compute-prior-dist`` or ``create-priors``); both are required together, exactly as FIMO
    requires. ``prior_alpha`` is the fraction of all binding sites that belong to the transcription factor of interest
    (0 to 1; FIMO ``--alpha``). The score of each occurrence is then a log-posterior odds score.

    ``p_value_threshold`` is the output threshold (FIMO ``--thresh``); with ``use_q_value_threshold`` it is applied to
    q-values instead (not combinable with ``compute_q_values=False``). ``scan_reverse_complement=False`` is
    ``--norc``; ``max_strand`` reports only the better strand where both pass. ``motif_ids`` selects motifs
    (``--motif``, repeated). ``background_file`` is a background model file, or give ``background_source`` as one of
    ``motif`` (``--bfile --motif--``), ``uniform`` or ``nrdb``. ``parse_genomic_coordinates=False`` is ``--no-pgc``.

    The new ``output_dir`` receives ``fimo.tsv``, ``fimo.gff``, ``fimo.html``, ``fimo.xml``, ``cisml.xml`` and
    ``best_site.narrowPeak``. Result keys: ``output_dir``, ``tsv_file``, ``gff_file``, ``match_count``, ``matches``
    (the first ``max_rows`` rows with ``motif_id``, ``motif_alt_id``, ``sequence_name``, ``start``, ``stop``,
    ``strand``, ``score``, ``p_value``, ``q_value``, ``matched_sequence``) and ``matches_truncated``.
    """
    motifs = _cli.existing_file(motif_file, "motif_file")
    sequences = _cli.existing_file(sequence_file, "sequence_file")
    prior = _cli.existing_file(prior_file, "prior_file")
    distribution = _cli.existing_file(prior_distribution_file, "prior_distribution_file")
    alpha = _cli.number(prior_alpha, "prior_alpha", minimum=0.0, maximum=1.0)
    threshold = _cli.number(p_value_threshold, "p_value_threshold", minimum=0.0, maximum=1.0)
    pseudocount = _cli.number(motif_pseudocount, "motif_pseudocount", minimum=0.0)
    stored = _cli.integer(max_stored_scores, "max_stored_scores", minimum=1)
    rows_limit = _cli.integer(max_rows, "max_rows", minimum=0)
    if not alpha > 0:
        raise ConfigurationError("prior_alpha must be greater than 0 and at most 1.")
    if use_q_value_threshold and not compute_q_values:
        raise ConfigurationError("use_q_value_threshold requires compute_q_values=True.")
    if background_file is not None and background_source is not None:
        raise ConfigurationError("Give either background_file or background_source, not both.")
    arguments: list[str] = [
        "--psp",
        str(prior),
        "--prior-dist",
        str(distribution),
        "--alpha",
        format(alpha, ".12g"),
    ]
    arguments += ["--thresh", format(threshold, ".12g")]
    arguments += ["--motif-pseudo", format(pseudocount, ".12g"), "--max-stored-scores", str(stored)]
    if use_q_value_threshold:
        arguments.append("--qv-thresh")
    if not compute_q_values:
        arguments.append("--no-qvalue")
    if not scan_reverse_complement:
        arguments.append("--norc")
    if max_strand:
        arguments.append("--max-strand")
    if not parse_genomic_coordinates:
        arguments.append("--no-pgc")
    if background_source is not None:
        arguments += [
            "--bfile",
            _BACKGROUND_SOURCES[
                _cli.choice(background_source, tuple(_BACKGROUND_SOURCES), "background_source")
            ],
        ]
    elif background_file is not None:
        arguments += ["--bfile", str(_cli.existing_file(background_file, "background_file"))]
    arguments += _motif_ids(motif_ids)
    out = Path(output_dir).expanduser().resolve()
    if out.exists():
        raise ConfigurationError(f"output_dir must not exist yet: {out}")
    executable = _cli.program("fimo", explicit=fimo_path, hint=_HINT)
    out.mkdir(parents=True)
    arguments += ["--verbosity", "1", "--oc", str(out), str(motifs), str(sequences)]
    expected = [out / "fimo.tsv", out / "fimo.gff"]
    result = _cli.run(
        "fimo", arguments, cwd=out, expected=expected, explicit=executable, hint=_HINT
    )
    result["output_dir"] = str(out)
    result["tsv_file"] = str(out / "fimo.tsv") if (out / "fimo.tsv").is_file() else None
    result["gff_file"] = str(out / "fimo.gff") if (out / "fimo.gff").is_file() else None
    result["outputs"] = [str(out / name) for name in _OUTPUT_FILES if (out / name).is_file()]
    if result["tsv_file"] is not None:
        total, rows = _parse_tsv(out / "fimo.tsv", rows_limit)
        result["match_count"] = total
        result["matches"] = rows
        result["matches_truncated"] = total > len(rows)
    else:
        result["match_count"] = 0
        result["matches"] = []
        result["matches_truncated"] = False
    return result


__all__ = ["scan_motifs_with_position_priors"]
_register_public_verbose(globals())
