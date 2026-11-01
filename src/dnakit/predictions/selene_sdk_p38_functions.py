"""selene-sdk sampling, feature look-up, genome window, ISM and strand-symmetric prediction operations (pack P38).

* ``sample_selene_dna_windows``: labelled DNA windows from random genome positions, from user-defined intervals,
  from a BED file or from a MAT/HDF5 matrix file, optionally exporting the sampled positions to a BED file.
* ``lookup_selene_genomic_features``: feature labels overlapping an interval, and the positive-feature test.
* ``get_selene_genome_window``: DNA sequence or one-hot encoding of a genomic interval.
* ``read_selene_ism_scores``: load an in silico mutagenesis result table and extract a per-base score matrix.
* ``predict_selene_sequences_strand_symmetric``: forward plus reverse-complement prediction with a torch model.

Install with ``python -m pip install "dnakit[bio-selene]"``. Feature BED files must be bgzip-compressed and tabix
indexed (``<name>.bed.gz`` and ``<name>.bed.gz.tbi``), as selene requires. No network is used. Tested with
selene-sdk 0.6.0.
"""

from __future__ import annotations

import os
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics import _p38_helpers as h
from dnakit.bioinformatics._common import require
from dnakit.exceptions import BackendExecutionError, ConfigurationError

_EXTRA = "bio-selene"
_SAMPLERS = ("random_positions", "intervals", "bed_file", "mat_file")
_MODES = ("train", "validate", "test")
_MAX_BATCH = 2048
_MAX_SHOWN = 50
_MAX_WINDOW = 100_000
_MAX_SEQS = 256
_BASES = "ACGT"


def _np() -> Any:
    return require("numpy", _EXTRA)


def _genome(path: Path) -> Any:
    sequences = require("selene_sdk.sequences", _EXTRA)
    try:
        return sequences.Genome(str(path))
    except (OSError, ValueError, KeyError) as exc:
        raise ConfigurationError(f"selene could not open the reference FASTA: {exc}") from exc


def _decode(window: Any) -> str:
    np = _np()
    out = []
    for row in window:
        out.append("N" if float(np.sum(row)) < 0.99 or np.max(row) < 0.99 else _BASES[int(np.argmax(row))])
    return "".join(out)


def _chrom_list(value: Sequence[str] | None, name: str) -> list[str]:
    if value is None:
        return []
    if isinstance(value, (str, bytes)):
        raise ConfigurationError(f"{name} must be a list of chromosome names.")
    return [h.text(c, name) for c in value]


def _feature_names(value: Sequence[str]) -> list[str]:
    if isinstance(value, (str, bytes)) or not value:
        raise ConfigurationError("feature_names must be a non-empty list of feature names.")
    return [h.text(f, "feature name") for f in value]


def _features(path: Path, names: list[str], threshold: Any) -> Any:
    targets = require("selene_sdk.targets", _EXTRA)
    if not Path(str(path) + ".tbi").is_file():
        raise ConfigurationError(f"The tabix index {path}.tbi is missing; bgzip and tabix-index the feature BED.")
    try:
        return targets.GenomicFeatures(str(path), names, feature_thresholds=threshold)
    except Exception as exc:  # selene raises tabix errors and ValueError
        raise ConfigurationError(f"selene could not open the feature file: {exc}") from exc


def sample_selene_dna_windows(
    sampler: str,
    *,
    reference_fasta: str | None = None,
    features_bed: str | None = None,
    feature_names: Sequence[str] | None = None,
    intervals_path: str | None = None,
    bed_path: str | None = None,
    bed_has_targets: bool = False,
    mat_path: str | None = None,
    sequence_key: str | None = None,
    targets_key: str | None = None,
    sequence_length: int = 1000,
    center_bin_to_predict: int = 200,
    feature_threshold: float = 0.5,
    batch_size: int = 8,
    mode: str = "train",
    validation_holdout: Sequence[str] = ("chr6", "chr7"),
    test_holdout: Sequence[str] = ("chr8", "chr9"),
    seed: int = 436,
    export_dir: str | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Draw labelled DNA windows with a selene sampler and optionally export the sampled positions.

    ``sampler`` selects the selene class: ``random_positions`` (``RandomPositionsSampler``: random positions
    of the genome), ``intervals`` (``IntervalsSampler``: positions inside the BED ``intervals_path``),
    ``bed_file`` (``BedFileSampler``: windows listed in ``bed_path``; with ``bed_has_targets`` the last column
    holds ``;``-separated feature indices) or ``mat_file`` (``MatFileSampler``: sequences and targets from a
    MAT or HDF5 file, keys ``sequence_key`` and ``targets_key``).

    The first two need ``reference_fasta``, a tabix-indexed ``features_bed`` (bgzip) and ``feature_names``;
    chromosomes in ``validation_holdout`` and ``test_holdout`` are kept out of ``train`` mode (``mode`` is
    ``train``, ``validate`` or ``test``). ``feature_threshold`` is the overlap fraction that makes a label
    positive. With ``export_dir`` the sampled positions are written to ``<export_dir>/<mode>_data.bed``
    (``OnlineSampler.save_dataset_to_file``; ``random_positions`` and ``intervals`` only). Returns the shapes,
    the decoded windows (at most 50) and their labels.
    """

    kind = h.choice(sampler, "sampler", _SAMPLERS)
    split = h.choice(mode, "mode", _MODES)
    batch = h.integer(batch_size, "batch_size", 1, _MAX_BATCH)
    length = h.integer(sequence_length, "sequence_length", 2, _MAX_WINDOW)
    center = h.integer(center_bin_to_predict, "center_bin_to_predict", 1, length)
    threshold = h.number(feature_threshold, "feature_threshold", 0.0, 1.0)
    seed_value = h.integer(seed, "seed", 0, 2**31 - 1)
    samplers = require("selene_sdk.samplers", _EXTRA)
    file_samplers = require("selene_sdk.samplers.file_samplers", _EXTRA)
    np = _np()
    handles: list[Any] = []
    export_target = None
    try:
        if kind in ("random_positions", "intervals"):
            if reference_fasta is None or features_bed is None or feature_names is None:
                raise ConfigurationError(f"{kind} needs reference_fasta, features_bed and feature_names.")
            fasta = h.existing(reference_fasta, "reference_fasta")
            names = _feature_names(feature_names)
            genome = _genome(fasta)
            targets = _features(h.existing(features_bed, "features_bed"), names, threshold)
            if export_dir is not None:
                export_target = Path(export_dir).expanduser()
                if not export_target.is_dir():
                    raise ConfigurationError("export_dir must be an existing directory.")
                if (export_target / f"{split}_data.bed").exists():
                    raise ConfigurationError(f"{split}_data.bed already exists in export_dir; refusing to overwrite.")
            common = dict(
                reference_sequence=genome, target_path=targets, features=names, seed=seed_value,
                validation_holdout=_chrom_list(validation_holdout, "validation_holdout"),
                test_holdout=_chrom_list(test_holdout, "test_holdout"),
                sequence_length=length, center_bin_to_predict=center, feature_thresholds=threshold, mode=split,
                save_datasets=[split] if export_target is not None else [],
                output_dir=str(export_target) if export_target is not None else None,
            )
            if kind == "random_positions":
                instance = samplers.RandomPositionsSampler(**common)
            else:
                if intervals_path is None:
                    raise ConfigurationError("intervals needs intervals_path.")
                instance = samplers.IntervalsSampler(intervals_path=str(h.existing(intervals_path, "intervals_path")), **common)
            sequences, labels = instance.sample(batch_size=batch)
            if export_target is not None:
                instance.save_dataset_to_file(split, close_filehandle=True)
        elif kind == "bed_file":
            if reference_fasta is None or bed_path is None:
                raise ConfigurationError("bed_file needs reference_fasta and bed_path.")
            if export_dir is not None:
                raise ConfigurationError("export_dir applies to random_positions and intervals only.")
            genome = _genome(h.existing(reference_fasta, "reference_fasta"))
            source = h.existing(bed_path, "bed_path")
            with source.open() as handle:
                rows = sum(1 for _ in handle)
            if rows == 0:
                raise ConfigurationError("bed_path is empty.")
            n_features = None if feature_names is None else len(_feature_names(feature_names))
            if bed_has_targets and not n_features:
                raise ConfigurationError("bed_has_targets needs feature_names to give the number of features.")
            instance = file_samplers.BedFileSampler(
                str(source), genome, rows, sequence_length=length, targets_avail=h.flag(bed_has_targets, "bed_has_targets"),
                n_features=n_features,
            )
            handles.append(instance)
            result = instance.sample(batch_size=batch)
            sequences, labels = (result[0], result[1]) if isinstance(result, tuple) and len(result) > 1 else (result[0], None)
        else:
            if mat_path is None or sequence_key is None:
                raise ConfigurationError("mat_file needs mat_path and sequence_key.")
            if export_dir is not None:
                raise ConfigurationError("export_dir applies to random_positions and intervals only.")
            instance = file_samplers.MatFileSampler(
                str(h.existing(mat_path, "mat_path")), h.text(sequence_key, "sequence_key"),
                targets_key=None if targets_key is None else h.text(targets_key, "targets_key"), random_seed=seed_value,
            )
            result = instance.sample(batch_size=batch)
            sequences = result[0]
            labels = result[1] if len(result) > 1 else None
    except ConfigurationError:
        raise
    except (KeyError, ValueError, TypeError, OSError, AttributeError, IndexError) as exc:
        raise ConfigurationError(f"selene could not sample windows: {exc}") from exc
    finally:
        for item in handles:
            handle = getattr(item, "_file_handle", None)
            if handle is not None:
                handle.close()
    arr = np.asarray(sequences)
    shown = min(len(arr), _MAX_SHOWN)
    exported = None
    if export_target is not None:
        exported = str(export_target / f"{split}_data.bed")
    return {
        "status": "completed",
        "operation": "sample_selene_dna_windows",
        "sampler": kind,
        "mode": split,
        "sequences_shape": [int(v) for v in arr.shape],
        "targets_shape": None if labels is None else [int(v) for v in np.asarray(labels).shape],
        "windows": [_decode(arr[i]) if arr.ndim == 3 and arr.shape[2] == 4 else None for i in range(shown)],
        "targets": None if labels is None else [[float(v) for v in row] for row in np.asarray(labels)[:shown]],
        "exported_bed": exported,
        "exported_bed_exists": bool(exported and os.path.isfile(exported)),
        "provenance": h.provenance("selene-sdk", f"selene_sdk.samplers ({kind})"),
    }


def lookup_selene_genomic_features(
    features_bed: str,
    feature_names: Sequence[str],
    chrom: str,
    start: int,
    end: int,
    *,
    strand: str | None = None,
    feature_threshold: float | None = 0.5,
    verbose: bool = True,
) -> dict[str, Any]:
    """Look up which genomic features overlap an interval, and whether any feature is positive there.

    Uses ``selene_sdk.targets.GenomicFeatures.get_feature_data`` (a 0/1 label per feature name for the interval
    ``chrom:start-end``, 0-based half-open) and ``is_positive``. ``features_bed`` is a bgzip-compressed,
    tabix-indexed BED with the feature name in column 4. A feature is labelled 1 when it overlaps at least
    ``feature_threshold`` of the interval (a fraction in [0, 1]; null requires any overlap), and the interval is
    positive when at least one feature is. ``strand`` is ``+``, ``-`` or null.
    """

    path = h.existing(features_bed, "features_bed")
    names = _feature_names(feature_names)
    name = h.text(chrom, "chrom")
    lo = h.integer(start, "start", 0, 4_000_000_000)
    hi = h.integer(end, "end", lo + 1, 4_000_000_000)
    side = None if strand is None else h.choice(strand, "strand", ("+", "-"))
    threshold = None if feature_threshold is None else h.number(feature_threshold, "feature_threshold", 0.0, 1.0)
    use = 0.0 if threshold is None else threshold
    targets = _features(path, names, use)
    try:
        data = targets.get_feature_data(name, lo, hi, strand=side)
        positive = bool(targets.is_positive(name, lo, hi, strand=side))
    except Exception as exc:  # tabix and selene raise several unrelated types
        raise ConfigurationError(f"selene could not query the features: {exc}") from exc
    labels = [int(v) for v in data]
    return {
        "status": "completed",
        "operation": "lookup_selene_genomic_features",
        "interval": {"chrom": name, "start": lo, "end": hi, "strand": side},
        "feature_threshold": threshold,
        "features": dict(zip(names, labels)),
        "overlapping_features": [n for n, v in zip(names, labels) if v],
        "is_positive": positive,
        "provenance": h.provenance("selene-sdk", "selene_sdk.targets.GenomicFeatures.get_feature_data / is_positive"),
    }


def get_selene_genome_window(
    reference_fasta: str,
    chrom: str,
    start: int,
    end: int,
    *,
    output: str = "sequence",
    strand: str = "+",
    pad: bool = False,
    verbose: bool = True,
) -> dict[str, Any]:
    """Get the DNA sequence or the one-hot encoding of a genomic interval from a FASTA, using
    ``selene_sdk.sequences.Genome``.

    ``output`` is ``sequence`` (``get_sequence_from_coords``: the upper-case sequence; the reverse complement when
    ``strand`` is ``-``) or ``encoding`` (``get_encoding_from_coords_check_unk``: a length-by-4 one-hot matrix in
    the order A, C, G, T, with an all-zero row per N, plus a flag telling whether unknown bases occurred).
    Coordinates are 0-based, half-open. With ``pad`` out-of-range bases become N; otherwise an out-of-range
    interval returns an empty result.
    """

    fasta = h.existing(reference_fasta, "reference_fasta")
    mode = h.choice(output, "output", ("sequence", "encoding"))
    name = h.text(chrom, "chrom")
    lo = h.integer(start, "start", 0, 4_000_000_000)
    hi = h.integer(end, "end", lo + 1, lo + _MAX_WINDOW)
    side = h.choice(strand, "strand", ("+", "-", "."))
    padded = h.flag(pad, "pad")
    genome = _genome(fasta)
    try:
        if mode == "sequence":
            text = genome.get_sequence_from_coords(name, lo, hi, strand=side, pad=padded)
            payload: dict[str, Any] = {"sequence": text, "length": len(text)}
        else:
            matrix, unknown = genome.get_encoding_from_coords_check_unk(name, lo, hi, strand=side, pad=padded)
            rows = [[float(v) for v in row] for row in matrix]
            payload = {"encoding": rows, "length": len(rows), "contains_unknown_bases": bool(unknown),
                       "alphabet_order": list(_BASES)}
    except ValueError as exc:
        raise ConfigurationError(f"selene could not read the interval: {exc}") from exc
    return {
        "status": "completed",
        "operation": "get_selene_genome_window",
        "interval": {"chrom": name, "start": lo, "end": hi, "strand": side},
        "output": mode,
        **payload,
        "provenance": h.provenance("selene-sdk", "selene_sdk.sequences.Genome"),
    }


def read_selene_ism_scores(
    ism_path: str,
    *,
    feature: str | None = None,
    reference_mask: float | None = None,
    verbose: bool = True,
) -> dict[str, Any]:
    """Load an in silico mutagenesis result table and extract a per-base score matrix for one feature.

    Uses ``selene_sdk.interpret.ISMResult.from_file`` (a tab-separated table with columns ``pos``, ``ref``,
    ``alt`` and one column per feature; an optional first row with ``ref`` and ``alt`` equal to ``NA`` holds the
    reference prediction) and ``ISMResult.get_score_matrix_for``. Returns the reference sequence and the feature
    names; with ``feature`` also the length-by-4 matrix (A, C, G, T) where alternative-base cells hold the
    mutant score and reference-base cells hold the reference prediction (or ``reference_mask``).
    """

    path = h.existing(ism_path, "ism_path")
    interpret = require("selene_sdk.interpret", _EXTRA)
    pandas = require("pandas", _EXTRA)
    try:
        result = interpret.ISMResult.from_file(str(path))
        columns = [str(c) for c in pandas.read_csv(path, sep="\t", nrows=0).columns if c not in ("pos", "ref", "alt")]
    except (ValueError, KeyError, OSError, Exception) as exc:  # selene raises bare Exception for incomplete tables
        raise ConfigurationError(f"selene could not read the ISM table: {exc}") from exc
    out: dict[str, Any] = {
        "status": "completed",
        "operation": "read_selene_ism_scores",
        "reference_sequence": result.reference_sequence,
        "length": len(result.reference_sequence),
        "features": columns,
        "alphabet_order": list(_BASES),
        "provenance": h.provenance("selene-sdk", "selene_sdk.interpret.ISMResult"),
    }
    if feature is not None:
        name = h.choice(feature, "feature", tuple(columns))
        mask = None if reference_mask is None else h.number(reference_mask, "reference_mask", -1e12, 1e12)
        try:
            matrix = result.get_score_matrix_for(name, reference_mask=mask)
        except (ValueError, KeyError) as exc:
            raise ConfigurationError(f"selene could not extract the score matrix: {exc}") from exc
        out["feature"] = name
        out["score_matrix"] = [[h.plain(float(v)) for v in row] for row in matrix]
    return out


def predict_selene_sequences_strand_symmetric(
    model_path: str,
    sequences: Sequence[str],
    *,
    mode: str = "mean",
    verbose: bool = True,
) -> dict[str, Any]:
    """Strand-symmetric prediction: combine a model's forward and reverse-complement outputs, using
    ``selene_sdk.utils.NonStrandSpecific``.

    ``model_path`` holds a whole ``torch.nn.Module`` saved with ``torch.save(model)`` (it is unpickled, so only
    use trusted files). Each sequence is one-hot encoded, the model is run on it and on its reverse complement,
    and the two outputs are combined with ``mode`` ``mean`` or ``max``. Sequences share one length. Returns the
    combined predictions and the plain forward predictions for comparison.
    """

    path = h.existing(model_path, "model_path")
    how = h.choice(mode, "mode", ("mean", "max"))
    if isinstance(sequences, (str, bytes)) or not sequences or len(sequences) > _MAX_SEQS:
        raise ConfigurationError(f"sequences must be a list of 1 to {_MAX_SEQS} DNA strings.")
    items = []
    for item in sequences:
        text = h.text(item, "sequence").upper()
        if set(text) - set("ACGTN"):
            raise ConfigurationError("sequences may contain only A, C, G, T or N.")
        items.append(text)
    if len({len(s) for s in items}) != 1 or len(items[0]) > _MAX_WINDOW:
        raise ConfigurationError("sequences must have the same length (at most 100000).")
    torch = require("torch", _EXTRA)
    np = _np()
    utils = require("selene_sdk.utils", _EXTRA)
    genome = require("selene_sdk.sequences", _EXTRA).Genome
    predict_module = require("selene_sdk.predict.model_predict", _EXTRA)
    try:
        model = torch.load(path, map_location="cpu", weights_only=False)
    except Exception as exc:  # torch raises several unrelated types for unreadable files
        raise ConfigurationError(f"model file could not be loaded: {exc}") from exc
    if not isinstance(model, torch.nn.Module):
        raise ConfigurationError("model file must hold a whole torch.nn.Module (torch.save(model)).")
    model.eval()
    batch = np.stack([np.asarray(genome.sequence_to_encoding(s), dtype=np.float32) for s in items])
    try:
        symmetric = utils.NonStrandSpecific(model, mode=how)
        symmetric.eval()
        combined = np.asarray(predict_module.predict(symmetric, batch, use_cuda=False))
        forward = np.asarray(predict_module.predict(model, batch, use_cuda=False))
    except (RuntimeError, ValueError, TypeError) as exc:
        raise BackendExecutionError(f"selene prediction failed: {exc}") from exc
    return {
        "status": "completed",
        "operation": "predict_selene_sequences_strand_symmetric",
        "mode": how,
        "sequence_count": len(items),
        "length": len(items[0]),
        "shape": [int(v) for v in combined.shape],
        "predictions": [[h.plain(float(v)) for v in row] for row in combined.reshape(len(items), -1)],
        "forward_only_predictions": [[h.plain(float(v)) for v in row] for row in forward.reshape(len(items), -1)],
        "provenance": h.provenance("selene-sdk", "selene_sdk.utils.NonStrandSpecific", device="cpu"),
    }


__all__ = [
    "get_selene_genome_window",
    "lookup_selene_genomic_features",
    "predict_selene_sequences_strand_symmetric",
    "read_selene_ism_scores",
    "sample_selene_dna_windows",
]

_register_public_verbose(globals(), tuple(__all__))
