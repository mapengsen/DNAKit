"""Optional regulatory-genomics adapters for motif scans, PyTorch models and fragment matrices.

The functions in this module are explicit adapters over optional packages:

- ``read_motif_file()``, ``scan_motif_hits()`` and ``score_motif_kmers()`` use
  ``gimmemotifs`` motif objects. They scan only the sequences passed in; no
  genome, background file or motif database is downloaded.
- ``predict_sequences()``, ``marginalize_motif_effect()`` and
  ``deep_lift_shap_attribution()`` call ``tangermeme`` on a PyTorch model that
  the caller supplies in eval mode. Models are never trained or loaded here.
- ``read_fragment_qc()``, ``count_tile_matrix()`` and ``count_peak_fragments()``
  import a five-column paired-end fragment file with ``snapatac2``.

Packages are imported only when a function runs. A missing package raises
``BackendUnavailableError`` with the ``pip install "dnakit[...]"`` hint.
"""

from __future__ import annotations

import importlib.metadata
import math
import numbers
import tempfile
from collections.abc import Iterator, Mapping, Sequence
from contextlib import contextmanager
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import input_path, named_sequence_items, require, sequence_text
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_GIMMEMOTIFS_EXTRA = "bio-gimmemotifs"
_TANGERMEME_EXTRA = "bio-tangermeme"
_SNAPATAC2_EXTRA = "bio-snapatac2"
_BASES = "ACGT"
_MOTIF_FORMATS = ("pfm", "transfac", "jaspar")
_COUNTING_STRATEGIES = ("fragment", "insertion", "paired-insertion")
_DEFAULT_EXCLUDED_CHROMS = ("chrM", "chrY", "M", "Y")
_DEFAULT_MAX_SEQUENCES = 1_000
_DEFAULT_MAX_SEQUENCE_LENGTH = 100_000
_MAX_SEQUENCES_LIMIT = 100_000
_MAX_SEQUENCE_LENGTH_LIMIT = 10_000_000
_MAX_MOTIFS = 1_000
_MAX_MOTIF_LENGTH = 100
_DEFAULT_MAX_MOTIF_BYTES = 10_000_000
_MAX_HITS_LIMIT = 100_000
_MAX_KMERS = 100_000
_DEFAULT_MAX_FRAGMENTS = 5_000_000
_DEFAULT_MAX_CELLS = 100_000
_DEFAULT_MAX_ENTRIES = 1_000_000
_MAX_ENTRIES_LIMIT = 50_000_000
_DEFAULT_MAX_PEAKS = 100_000
_MAX_PEAKS_LIMIT = 10_000_000
_MAX_BIN_SIZE = 10_000_000
_DEFAULT_MAX_TENSOR_SEQUENCES = 256
_DEFAULT_MAX_TENSOR_LENGTH = 10_000
_DEFAULT_MAX_POSITIONS = 200_000
_MAX_POSITIONS_LIMIT = 5_000_000
_FRAGMENT_COLUMNS = 5
_FRAGMENT_QC_COLUMNS = ("n_fragment", "frac_dup", "frac_mito")


def _integer(value: object, name: str, low: int, high: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not low <= value <= high:
        raise ConfigurationError(f"{name} must be an integer in [{low}, {high}].")
    return value


def _fraction(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ConfigurationError(f"{name} must be a number in [0, 1].")
    number = float(value)
    if not math.isfinite(number) or not 0.0 <= number <= 1.0:
        raise ConfigurationError(f"{name} must be a number in [0, 1].")
    return number


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be True or False.")
    return value


def _label(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.strip() or "\x00" in value:
        raise ConfigurationError(f"{name} must be a non-empty string without NUL characters.")
    return value


def _choose(value: object, allowed: Sequence[str], name: str) -> str:
    if value not in allowed:
        raise ConfigurationError(f"{name} must be one of {', '.join(allowed)}; got {value!r}.")
    return str(value)


def _finite(value: object) -> float | None:
    number = float(value)  # type: ignore[arg-type]
    return number if math.isfinite(number) else None


def _nested(value: Any) -> Any:
    """Replace non-finite floats with None inside nested lists."""

    if isinstance(value, list):
        return [_nested(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _version(package: str) -> str:
    try:
        return importlib.metadata.version(package)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _dna(value: object, name: str, max_length: int) -> str:
    text = sequence_text(value)  # type: ignore[arg-type]
    if len(text) > max_length:
        raise ConfigurationError(f"{name} is longer than the limit of {max_length} bases.")
    if set(text) - set(_BASES):
        raise ConfigurationError(f"{name} contains non-ACGT symbols; use A, C, G and T only.")
    return text


def _named_dna(sequences: object, max_sequences: int, max_length: int) -> list[tuple[str, str]]:
    items = named_sequence_items(sequences)  # type: ignore[arg-type]
    if not items:
        raise ConfigurationError("At least one sequence is required.")
    if len(items) > max_sequences:
        raise ConfigurationError(
            f"sequence count exceeds max_sequences={max_sequences}; increase the explicit limit."
        )
    return [
        (_label(seq_id, "sequence ID"), _dna(value, f"sequence {seq_id!r}", max_length))
        for seq_id, value in items
    ]


def _pfm(value: object, name: str) -> list[list[float]]:
    if (
        isinstance(value, (str, bytes))
        or not isinstance(value, Sequence)
        or not 1 <= len(value) <= _MAX_MOTIF_LENGTH
    ):
        raise ConfigurationError(
            f"{name} must be a list of 1 to {_MAX_MOTIF_LENGTH} rows of A, C, G and T counts."
        )
    rows: list[list[float]] = []
    for row in value:
        if isinstance(row, (str, bytes)) or not isinstance(row, Sequence) or len(row) != 4:
            raise ConfigurationError(f"{name} rows must contain exactly four A, C, G, T counts.")
        counts: list[float] = []
        for count in row:
            if (
                isinstance(count, bool)
                or not isinstance(count, numbers.Real)
                or not math.isfinite(count)
                or count < 0
            ):
                raise ConfigurationError(f"{name} counts must be finite non-negative numbers.")
            counts.append(float(count))
        if sum(counts) <= 0:
            raise ConfigurationError(f"{name} rows must not be all zero.")
        rows.append(counts)
    return rows


def _motif_set(motifs: object, max_motifs: int) -> list[tuple[str, list[list[float]]]]:
    if not isinstance(motifs, Mapping) or not motifs:
        raise ConfigurationError("motifs must be a non-empty mapping of motif ID to PFM counts.")
    if len(motifs) > max_motifs:
        raise ConfigurationError(
            f"motif count exceeds max_motifs={max_motifs}; increase the explicit limit."
        )
    return [
        (_label(key, "motif ID"), _pfm(matrix, f"motif {key!r}")) for key, matrix in motifs.items()
    ]


def _upstream(call: Any) -> Any:
    """Run an upstream call and convert its failures into DNAKit errors."""

    try:
        return call()
    except (KeyboardInterrupt, SystemExit):
        raise
    except BaseException as exc:  # pyo3 panics derive from BaseException, not Exception
        raise BackendExecutionError(f"snapatac2 failed: {exc}") from exc


def _chrom_sizes(value: object) -> dict[str, int]:
    if not isinstance(value, Mapping) or not value:
        raise ConfigurationError("chrom_sizes must be a non-empty mapping of chromosome to length.")
    sizes: dict[str, int] = {}
    for name, length in value.items():
        chrom = _label(name, "chromosome name")
        sizes[chrom] = _integer(length, f"length of {chrom}", 1, 10**12)
    return sizes


def _tab_fields(line: str, source: str, number: int) -> list[str]:
    if line.startswith("#") or not line.strip():
        raise ConfigurationError(f"{source} line {number} is empty or a comment; remove it.")
    return line.rstrip("\r\n").split("\t")


def _span(fields: list[str], sizes: Mapping[str, int], source: str, number: int) -> None:
    chrom = fields[0]
    if chrom not in sizes:
        raise ConfigurationError(
            f"{source} line {number} uses chromosome {chrom!r}, not in chrom_sizes."
        )
    try:
        start, end = int(fields[1]), int(fields[2])
    except ValueError as exc:
        raise ConfigurationError(f"{source} line {number} has non-integer coordinates.") from exc
    if not 0 <= start < end <= sizes[chrom]:
        raise ConfigurationError(
            f"{source} line {number} must satisfy 0 <= start < end <= {sizes[chrom]} on {chrom}."
        )


def _fragment_file(value: object, sizes: Mapping[str, int], max_fragments: int) -> Path:
    """Validate a five-column paired-end fragment file without loading it into memory."""

    path = input_path(value)  # type: ignore[arg-type]
    count = 0
    with path.open("rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            fields = _tab_fields(line, path.name, number)
            if len(fields) != _FRAGMENT_COLUMNS:
                raise ConfigurationError(
                    f"{path.name} line {number} has {len(fields)} tab-separated columns; "
                    f"expected {_FRAGMENT_COLUMNS}."
                )
            _span(fields, sizes, path.name, number)
            if not fields[3]:
                raise ConfigurationError(f"{path.name} line {number} has an empty barcode.")
            try:
                reads = int(fields[4])
            except ValueError as exc:
                raise ConfigurationError(
                    f"{path.name} line {number} has a non-integer count."
                ) from exc
            if reads < 1:
                raise ConfigurationError(
                    f"{path.name} line {number} must have a count of at least 1."
                )
            count += 1
            if count > max_fragments:
                raise ConfigurationError(
                    f"{path.name} has more than max_fragments={max_fragments} rows; "
                    "increase the explicit limit."
                )
    if count == 0:
        raise ConfigurationError(f"{path.name} contains no fragments.")
    return path


def _peak_file(value: object, sizes: Mapping[str, int], max_peaks: int) -> Path:
    """Validate a BED file with at least three columns; only chrom, start and end are used."""

    path = input_path(value)  # type: ignore[arg-type]
    count = 0
    with path.open("rt", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            fields = _tab_fields(line, path.name, number)
            if len(fields) < 3:
                raise ConfigurationError(f"{path.name} line {number} needs chrom, start and end.")
            _span(fields, sizes, path.name, number)
            count += 1
            if count > max_peaks:
                raise ConfigurationError(
                    f"{path.name} has more than max_peaks={max_peaks} peaks; "
                    "increase the explicit limit."
                )
    if count == 0:
        raise ConfigurationError(f"{path.name} contains no peaks.")
    return path


@contextmanager
def _imported_fragments(
    snap: Any,
    fragment_file: object,
    chrom_sizes: object,
    *,
    min_num_fragments: int,
    max_fragments: int,
    max_cells: int,
) -> Iterator[tuple[Any, Path]]:
    """Import fragments into a temporary snapATAC2 file that is removed on exit."""

    sizes = _chrom_sizes(chrom_sizes)
    path = _fragment_file(fragment_file, sizes, max_fragments)
    with tempfile.TemporaryDirectory(prefix="dnakit-snapatac2-") as scratch:
        work = Path(scratch)
        adata = _upstream(
            lambda: snap.pp.import_fragments(
                path,
                sizes,
                is_paired=True,
                min_num_fragments=min_num_fragments,
                file=work / "fragments.h5ad",
                sorted_by_barcode=False,
                tempdir=work,
                n_jobs=1,
            )
        )
        try:
            if adata.n_obs > max_cells:
                raise ConfigurationError(
                    f"{path.name} yields {adata.n_obs} cells; max_cells={max_cells}."
                )
            yield adata, work
        finally:
            adata.close()


def _as_list(value: Any) -> list[Any]:
    if hasattr(value, "to_list"):
        return list(value.to_list())
    if hasattr(value, "tolist"):
        return list(value.tolist())
    return list(value)


def _count_entries(matrix: Any, max_entries: int) -> tuple[list[list[int | float]], int]:
    """Return sorted [row, column, value] triplets from a sparse or dense matrix."""

    np = require("numpy", _SNAPATAC2_EXTRA)
    if hasattr(matrix, "tocoo"):
        coo = matrix.tocoo()
        rows, cols, values = coo.row, coo.col, coo.data
    else:
        dense = np.asarray(matrix)
        rows, cols = np.nonzero(dense)
        values = dense[rows, cols]
    if len(values) > max_entries:
        raise ConfigurationError(
            f"the count matrix has {len(values)} non-zero entries; max_entries={max_entries}."
        )
    order = np.lexsort((cols, rows))
    triplets: list[list[int | float]] = []
    for index in order:
        value = float(values[index])
        triplets.append(
            [int(rows[index]), int(cols[index]), int(value) if value.is_integer() else value]
        )
    total = sum(item[2] for item in triplets)
    return triplets, int(total) if float(total).is_integer() else float(total)


def _snap_provenance(**notes: object) -> dict[str, Any]:
    return {
        "package": "snapatac2",
        "version": _version("snapatac2"),
        "input_format": "five tab-separated columns: chrom, start, end, barcode, read count",
        "coordinates": "0-based half-open fragment and bin coordinates, as named by snapatac2",
        **notes,
    }


def _gimme_provenance(**notes: object) -> dict[str, Any]:
    return {
        "package": "gimmemotifs",
        "version": _version("gimmemotifs"),
        **notes,
    }


def read_motif_file(
    path: str | Path,
    *,
    fmt: str = "pfm",
    max_bytes: int = _DEFAULT_MAX_MOTIF_BYTES,
    max_motifs: int = _MAX_MOTIFS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Parse PFM, TRANSFAC or JASPAR motifs with gimmemotifs.

    ``fmt`` is ``"pfm"``, ``"transfac"`` or ``"jaspar"``. Each motif is returned
    with its PFM counts, so the ``pfm`` value can be passed to ``scan_motif_hits``.
    """

    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    source = input_path(path)
    selected = _choose(fmt, _MOTIF_FORMATS, "fmt")
    size_limit = _integer(max_bytes, "max_bytes", 1, 100_000_000)
    motif_limit = _integer(max_motifs, "max_motifs", 1, _MAX_MOTIFS)
    size = source.stat().st_size
    if size > size_limit:
        raise ConfigurationError(f"{source.name} has {size} bytes; max_bytes={size_limit}.")
    try:
        parsed = motif_module.read_motifs(str(source), fmt=selected, as_dict=False)
    except (OSError, ValueError, TypeError, IndexError, KeyError) as exc:
        raise ConfigurationError(
            f"gimmemotifs could not read {source.name} as {selected}: {exc}"
        ) from exc
    if not parsed:
        raise ConfigurationError(f"{source.name} contains no {selected} motifs.")
    if len(parsed) > motif_limit:
        raise ConfigurationError(f"{source.name} has more than max_motifs={motif_limit} motifs.")
    motifs = []
    for motif in parsed:
        pfm = [[float(value) for value in row] for row in motif.pfm]
        motifs.append(
            {
                "id": str(motif.id),
                "length": len(pfm),
                "consensus": str(motif.consensus),
                "information_content": _finite(motif.information_content),
                "min_score": _finite(motif.min_score),
                "max_score": _finite(motif.max_score),
                "pfm": pfm,
            }
        )
    return {
        "status": "completed",
        "operation": "read_motif_file",
        "source": str(source),
        "format": selected,
        "motif_count": len(motifs),
        "motifs": motifs,
        "provenance": _gimme_provenance(
            counts="PFM counts as stored; gimmemotifs adds its pseudocount when scoring",
        ),
    }


def scan_motif_hits(
    sequences: object,
    motifs: Mapping[str, Sequence[Sequence[float]]],
    *,
    threshold: float = 0.9,
    scan_rc: bool = True,
    max_hits_per_sequence: int = 100,
    max_sequences: int = _DEFAULT_MAX_SEQUENCES,
    max_sequence_length: int = _DEFAULT_MAX_SEQUENCE_LENGTH,
    verbose: bool = True,
) -> dict[str, Any]:
    """Scan DNA sequences with PFM motifs using gimmemotifs log-odds scores.

    ``motifs`` maps motif IDs to PFM count rows. ``threshold`` is the cutoff as a
    fraction of each motif's own score range (minimum to maximum). The limit
    ``max_hits_per_sequence`` applies to each motif in each sequence.
    """

    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    named = _named_dna(
        sequences,
        _integer(max_sequences, "max_sequences", 1, _MAX_SEQUENCES_LIMIT),
        _integer(max_sequence_length, "max_sequence_length", 1, _MAX_SEQUENCE_LENGTH_LIMIT),
    )
    motif_items = _motif_set(motifs, _MAX_MOTIFS)
    cutoff = _fraction(threshold, "threshold")
    both_strands = _flag(scan_rc, "scan_rc")
    limit = _integer(max_hits_per_sequence, "max_hits_per_sequence", 1, _MAX_HITS_LIMIT)
    sequence_order = {seq_id: index for index, (seq_id, _text) in enumerate(named)}
    fasta = dict(named)
    hits: list[tuple[tuple[int, int, int, int], dict[str, Any]]] = []
    for motif_index, (motif_id, matrix) in enumerate(motif_items):
        motif = motif_module.Motif(matrix)
        found = motif.scan_all(fasta, cutoff=cutoff, nreport=limit + 1, scan_rc=both_strands)
        for seq_id, _text in named:
            matches = found.get(seq_id, [])
            if len(matches) > limit:
                raise ConfigurationError(
                    f"motif {motif_id!r} has more than max_hits_per_sequence={limit} hits in "
                    f"sequence {seq_id!r}; raise the limit or the threshold."
                )
            for start, score, strand in matches:
                sign = "+" if strand > 0 else "-"
                hit = {
                    "sequence": seq_id,
                    "motif": motif_id,
                    "start": int(start),
                    "end": int(start) + len(matrix),
                    "strand": sign,
                    "score": _finite(score),
                }
                key = (sequence_order[seq_id], motif_index, int(start), 0 if sign == "+" else 1)
                hits.append((key, hit))
    ordered = [hit for _key, hit in sorted(hits, key=lambda item: item[0])]
    return {
        "status": "completed",
        "operation": "scan_motif_hits",
        "motif_count": len(motif_items),
        "sequence_count": len(named),
        "threshold": cutoff,
        "hit_count": len(ordered),
        "hits": ordered,
        "provenance": _gimme_provenance(
            score="gimmemotifs log-odds: sum of ln(ppm / 0.25 + 0.01) over motif positions",
            cutoff="min + (max - min) * threshold, using each motif's own score range",
            coordinates=(
                "0-based start on the forward sequence; end is exclusive; "
                "minus-strand hits use the forward start"
            ),
        ),
    }


def score_motif_kmers(
    motif: Sequence[Sequence[float]],
    kmers: Sequence[str],
    *,
    max_kmers: int = 10_000,
    verbose: bool = True,
) -> dict[str, Any]:
    """Score equal-length DNA k-mers against one PFM with gimmemotifs log-odds."""

    motif_module = require("gimmemotifs.motif", _GIMMEMOTIFS_EXTRA)
    rows = _pfm(motif, "motif")
    limit = _integer(max_kmers, "max_kmers", 1, _MAX_KMERS)
    if isinstance(kmers, (str, bytes)) or not isinstance(kmers, Sequence) or not kmers:
        raise ConfigurationError("kmers must be a non-empty list of DNA strings.")
    if len(kmers) > limit:
        raise ConfigurationError(
            f"kmer count exceeds max_kmers={limit}; increase the explicit limit."
        )
    texts = []
    for index, value in enumerate(kmers):
        text = _dna(value, f"kmer {index}", len(rows))
        if len(text) != len(rows):
            raise ConfigurationError(
                f"kmer {text!r} must have length {len(rows)} to match the motif."
            )
        texts.append(text)
    model = motif_module.Motif(rows)
    scores = [{"kmer": text, "score": _finite(model.score_kmer(text))} for text in texts]
    return {
        "status": "completed",
        "operation": "score_motif_kmers",
        "motif_length": len(rows),
        "consensus": str(model.consensus),
        "kmer_count": len(scores),
        "max_score": _finite(model.max_score),
        "min_score": _finite(model.min_score),
        "scores": scores,
        "provenance": _gimme_provenance(
            score="gimmemotifs log-odds: sum of ln(ppm / 0.25 + 0.01) over motif positions",
        ),
    }


def read_fragment_qc(
    fragment_file: str | Path,
    chrom_sizes: Mapping[str, int],
    *,
    min_num_fragments: int = 1,
    max_fragments: int = _DEFAULT_MAX_FRAGMENTS,
    max_cells: int = _DEFAULT_MAX_CELLS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-barcode fragment QC (fragment count, duplicate and mitochondrial fractions).

    ``fragment_file`` is a plain-text paired-end file with five tab-separated
    columns: chromosome, start, end, barcode and read count. ``min_num_fragments``
    is passed to snapATAC2 and removes barcodes with fewer fragments.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    minimum = _integer(min_num_fragments, "min_num_fragments", 1, 1_000_000)
    with _imported_fragments(
        snap,
        fragment_file,
        chrom_sizes,
        min_num_fragments=minimum,
        max_fragments=_integer(max_fragments, "max_fragments", 1, 1_000_000_000),
        max_cells=_integer(max_cells, "max_cells", 1, 10_000_000),
    ) as (adata, _work):
        barcodes = [str(item) for item in adata.obs_names]
        columns = {name: _as_list(adata.obs[name][:]) for name in _FRAGMENT_QC_COLUMNS}
        cells = [
            {
                "barcode": barcode,
                "n_fragment": int(columns["n_fragment"][index]),
                "frac_dup": _finite(columns["frac_dup"][index]),
                "frac_mito": _finite(columns["frac_mito"][index]),
            }
            for index, barcode in enumerate(barcodes)
        ]
    return {
        "status": "completed",
        "operation": "read_fragment_qc",
        "cell_count": len(cells),
        "cells": cells,
        "provenance": _snap_provenance(
            n_fragment="fragments on non-mitochondrial chromosomes, as counted by snapatac2",
            frac_dup="duplicate reads divided by all reads (the count column)",
            frac_mito="fragments on chrM/M divided by all fragments",
        ),
    }


def count_tile_matrix(
    fragment_file: str | Path,
    chrom_sizes: Mapping[str, int],
    *,
    bin_size: int = 500,
    counting_strategy: str = "fragment",
    exclude_chroms: Sequence[str] = _DEFAULT_EXCLUDED_CHROMS,
    min_num_fragments: int = 1,
    max_fragments: int = _DEFAULT_MAX_FRAGMENTS,
    max_cells: int = _DEFAULT_MAX_CELLS,
    max_entries: int = _DEFAULT_MAX_ENTRIES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cell-by-bin count matrix from a paired-end fragment file, using snapATAC2 tiles.

    Bins are fixed-width windows named ``chrom:start-end``. ``counting_strategy``
    is ``"fragment"``, ``"insertion"`` or ``"paired-insertion"``. The matrix is
    returned as sorted ``[cell_index, bin_index, count]`` triplets.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    width = _integer(bin_size, "bin_size", 1, _MAX_BIN_SIZE)
    strategy = _choose(counting_strategy, _COUNTING_STRATEGIES, "counting_strategy")
    excluded = [_label(name, "exclude_chroms entry") for name in exclude_chroms]
    entry_limit = _integer(max_entries, "max_entries", 1, _MAX_ENTRIES_LIMIT)
    with _imported_fragments(
        snap,
        fragment_file,
        chrom_sizes,
        min_num_fragments=_integer(min_num_fragments, "min_num_fragments", 1, 1_000_000),
        max_fragments=_integer(max_fragments, "max_fragments", 1, 1_000_000_000),
        max_cells=_integer(max_cells, "max_cells", 1, 10_000_000),
    ) as (adata, _work):
        _upstream(
            lambda: snap.pp.add_tile_matrix(
                adata,
                bin_size=width,
                counting_strategy=strategy,
                exclude_chroms=excluded,
                n_jobs=1,
            )
        )
        barcodes = [str(item) for item in adata.obs_names]
        features = [str(item) for item in adata.var_names]
        entries, total = _count_entries(adata.X[:], entry_limit)
    return {
        "status": "completed",
        "operation": "count_tile_matrix",
        "bin_size": width,
        "counting_strategy": strategy,
        "shape": [len(barcodes), len(features)],
        "barcodes": barcodes,
        "features": features,
        "nonzero_count": len(entries),
        "total_count": total,
        "entries": entries,
        "provenance": _snap_provenance(
            tiles="fixed-width bins from snapatac2 add_tile_matrix",
            excluded_chromosomes=excluded,
        ),
    }


def count_peak_fragments(
    fragment_file: str | Path,
    chrom_sizes: Mapping[str, int],
    peak_file: str | Path,
    *,
    counting_strategy: str = "fragment",
    min_num_fragments: int = 1,
    max_fragments: int = _DEFAULT_MAX_FRAGMENTS,
    max_cells: int = _DEFAULT_MAX_CELLS,
    max_peaks: int = _DEFAULT_MAX_PEAKS,
    max_entries: int = _DEFAULT_MAX_ENTRIES,
    verbose: bool = True,
) -> dict[str, Any]:
    """Cell-by-peak count matrix from fragments and a BED file of peak regions.

    ``peak_file`` needs chromosome, start and end columns; extra columns are
    ignored. The matrix uses snapATAC2 ``make_peak_matrix`` and is returned as
    sorted ``[cell_index, peak_index, count]`` triplets.
    """

    snap = require("snapatac2", _SNAPATAC2_EXTRA)
    strategy = _choose(counting_strategy, _COUNTING_STRATEGIES, "counting_strategy")
    peak_limit = _integer(max_peaks, "max_peaks", 1, _MAX_PEAKS_LIMIT)
    entry_limit = _integer(max_entries, "max_entries", 1, _MAX_ENTRIES_LIMIT)
    sizes = _chrom_sizes(chrom_sizes)
    peaks = _peak_file(peak_file, sizes, peak_limit)
    with _imported_fragments(
        snap,
        fragment_file,
        sizes,
        min_num_fragments=_integer(min_num_fragments, "min_num_fragments", 1, 1_000_000),
        max_fragments=_integer(max_fragments, "max_fragments", 1, 1_000_000_000),
        max_cells=_integer(max_cells, "max_cells", 1, 10_000_000),
    ) as (adata, work):
        peak_matrix = _upstream(
            lambda: snap.pp.make_peak_matrix(
                adata,
                peak_file=peaks,
                inplace=False,
                file=work / "peaks.h5ad",
                counting_strategy=strategy,
                chunk_size=500,
            )
        )
        try:
            barcodes = [str(item) for item in peak_matrix.obs_names]
            features = [str(item) for item in peak_matrix.var_names]
            entries, total = _count_entries(peak_matrix.X[:], entry_limit)
        finally:
            peak_matrix.close()
    return {
        "status": "completed",
        "operation": "count_peak_fragments",
        "counting_strategy": strategy,
        "shape": [len(barcodes), len(features)],
        "barcodes": barcodes,
        "features": features,
        "nonzero_count": len(entries),
        "total_count": total,
        "entries": entries,
        "provenance": _snap_provenance(peaks=str(peaks.name)),
    }


def _torch_model(model: object) -> Any:
    torch = require("torch", _TANGERMEME_EXTRA)
    if not isinstance(model, torch.nn.Module):
        raise ConfigurationError("model must be a torch.nn.Module.")
    if model.training:
        raise ConfigurationError(
            "call model.eval() first; dropout and batch norm would make outputs depend on mode."
        )
    return torch


def _one_hot_batch(torch: Any, utils: Any, named: Sequence[tuple[str, str]]) -> Any:
    lengths = {len(text) for _seq_id, text in named}
    if len(lengths) != 1:
        raise ConfigurationError("sequences must all have the same length for a model batch.")
    return torch.stack([utils.one_hot_encode(text) for _seq_id, text in named]).float()


def _single_output(output: object) -> Any:
    if isinstance(output, (list, tuple)):
        raise ConfigurationError("models with several output tensors are not supported.")
    return output.detach().cpu()  # type: ignore[attr-defined]


def _tensor_provenance(torch: Any, **notes: object) -> dict[str, Any]:
    return {
        "package": "tangermeme",
        "version": _version("tangermeme"),
        "torch": torch.__version__,
        **notes,
    }


def predict_sequences(
    model: Any,
    sequences: object,
    *,
    batch_size: int = 32,
    max_sequences: int = _DEFAULT_MAX_TENSOR_SEQUENCES,
    max_sequence_length: int = _DEFAULT_MAX_TENSOR_LENGTH,
    verbose: bool = True,
) -> dict[str, Any]:
    """Run a single-output PyTorch model over one-hot DNA with tangermeme.predict.

    The model must already be in eval mode. Sequences must share one length and
    contain only A, C, G and T.
    """

    torch = _torch_model(model)
    predict_module = require("tangermeme.predict", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    named = _named_dna(
        sequences,
        _integer(max_sequences, "max_sequences", 1, _MAX_SEQUENCES_LIMIT),
        _integer(max_sequence_length, "max_sequence_length", 1, _MAX_SEQUENCE_LENGTH_LIMIT),
    )
    batch = _integer(batch_size, "batch_size", 1, 4096)
    inputs = _one_hot_batch(torch, utils, named)
    output = _single_output(
        _upstream(lambda: predict_module.predict(model, inputs, batch_size=batch))
    )
    values = output.tolist()
    return {
        "status": "completed",
        "operation": "predict_sequences",
        "sequence_count": len(named),
        "shape": list(output.shape),
        "predictions": _nested(values),
        "provenance": _tensor_provenance(torch, function="tangermeme.predict.predict"),
    }


def marginalize_motif_effect(
    model: Any,
    sequences: object,
    motif: str,
    *,
    start: int | None = None,
    max_sequences: int = _DEFAULT_MAX_TENSOR_SEQUENCES,
    max_sequence_length: int = _DEFAULT_MAX_TENSOR_LENGTH,
    verbose: bool = True,
) -> dict[str, Any]:
    """In-silico marginalization: model output before and after substituting a motif.

    The same ``motif`` (an ACGT string) is written over each sequence. When
    ``start`` is None, tangermeme places the motif at the centre of each sequence.
    ``delta`` is ``y_after - y_before`` for each sequence.
    """

    torch = _torch_model(model)
    marginalize_module = require("tangermeme.marginalize", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    named = _named_dna(
        sequences,
        _integer(max_sequences, "max_sequences", 1, _MAX_SEQUENCES_LIMIT),
        _integer(max_sequence_length, "max_sequence_length", 1, _MAX_SEQUENCE_LENGTH_LIMIT),
    )
    motif_text = _dna(motif, "motif", _MAX_MOTIF_LENGTH)
    length = len(named[0][1])
    if any(len(text) != length for _seq_id, text in named):
        raise ConfigurationError("sequences must all have the same length for a model batch.")
    if len(motif_text) > length:
        raise ConfigurationError("motif is longer than the sequences.")
    if start is not None:
        start = _integer(start, "start", 0, length - len(motif_text))
    used_start = start if start is not None else length // 2 - len(motif_text) // 2
    inputs = _one_hot_batch(torch, utils, named)
    motif_tensor = utils.one_hot_encode(motif_text).float()[None]
    result = _upstream(
        lambda: marginalize_module.marginalize(model, inputs, motif_tensor, start=start)
    )
    before = _single_output(result.y_before)
    after = _single_output(result.y_after)
    return {
        "status": "completed",
        "operation": "marginalize_motif_effect",
        "sequence_count": len(named),
        "motif": motif_text,
        "start": used_start,
        "shape": list(before.shape),
        "y_before": _nested(before.tolist()),
        "y_after": _nested(after.tolist()),
        "delta": _nested((after - before).tolist()),
        "provenance": _tensor_provenance(torch, function="tangermeme.marginalize.marginalize"),
    }


def deep_lift_shap_attribution(
    model: Any,
    sequences: object,
    *,
    target: int = 0,
    n_shuffles: int = 20,
    random_state: int = 0,
    max_sequences: int = _DEFAULT_MAX_TENSOR_SEQUENCES,
    max_sequence_length: int = _DEFAULT_MAX_TENSOR_LENGTH,
    max_positions: int = _DEFAULT_MAX_POSITIONS,
    verbose: bool = True,
) -> dict[str, Any]:
    """Per-base DeepLIFT/SHAP attributions for one model output, using tangermeme.

    References are dinucleotide-shuffled copies of each sequence, drawn with
    ``random_state``. ``attributions`` is indexed as sequence, position, base
    (A, C, G, T); the limit ``max_positions`` applies to sequences times length.
    """

    torch = _torch_model(model)
    dl_module = require("tangermeme.deep_lift_shap", _TANGERMEME_EXTRA)
    utils = require("tangermeme.utils", _TANGERMEME_EXTRA)
    named = _named_dna(
        sequences,
        _integer(max_sequences, "max_sequences", 1, _MAX_SEQUENCES_LIMIT),
        _integer(max_sequence_length, "max_sequence_length", 1, _MAX_SEQUENCE_LENGTH_LIMIT),
    )
    limit = _integer(max_positions, "max_positions", 1, _MAX_POSITIONS_LIMIT)
    target_index = _integer(target, "target", 0, 100_000)
    shuffles = _integer(n_shuffles, "n_shuffles", 1, 1_000)
    seed = _integer(random_state, "random_state", 0, 2**32 - 1)
    positions = len(named) * len(named[0][1])
    if positions > limit:
        raise ConfigurationError(f"sequences contain {positions} positions; max_positions={limit}.")
    inputs = _one_hot_batch(torch, utils, named)
    attributions = _upstream(
        lambda: dl_module.deep_lift_shap(
            model, inputs, target=target_index, n_shuffles=shuffles, random_state=seed
        )
    )
    values = _single_output(attributions).permute(0, 2, 1)
    return {
        "status": "completed",
        "operation": "deep_lift_shap_attribution",
        "sequence_count": len(named),
        "shape": list(values.shape),
        "axis_order": "sequence, position, base (A, C, G, T)",
        "attributions": _nested(values.tolist()),
        "provenance": _tensor_provenance(
            torch,
            function="tangermeme.deep_lift_shap.deep_lift_shap",
            references="dinucleotide-shuffled copies of each sequence",
            random_state=seed,
        ),
    }


__all__ = [
    "count_peak_fragments",
    "count_tile_matrix",
    "deep_lift_shap_attribution",
    "marginalize_motif_effect",
    "predict_sequences",
    "read_fragment_qc",
    "read_motif_file",
    "scan_motif_hits",
    "score_motif_kmers",
]

_register_public_verbose(globals(), tuple(__all__))
