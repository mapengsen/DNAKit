"""Offline AlphaGenome variant-score utilities (review group G18q).

Each public function imports ``alphagenome``, ``anndata`` and ``pandas`` only when it runs.
Inputs are local: serialized ``DenseVariantScores`` protobuf messages, AnnData score tables,
and track-metadata records. Model inference, Atlas and API clients, gRPC calls, and plotting
are not wrapped here.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from dnakit._verbosity import register_public as _register_public_verbose
from dnakit.bioinformatics._common import require
from dnakit.exceptions import ConfigurationError
from dnakit.predictions.dna_model_backends import (
    _choice,
    _clean,
    _provenance,
    _rejected,
    _text,
)

_ALPHAGENOME = "alphagenome"
_EXTRAS = "bio-alphagenome"
_STRANDS = ("+", "-", ".")
_ORGANISMS = ("homo_sapiens", "mus_musculus")
_MAX_ENTRIES = 10_000


def _flag(value: object, name: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{name} must be true or false.")
    return value


def _records_list(value: object, name: str) -> list[dict[str, Any]]:
    """A non-empty list of mappings, copied with string keys."""

    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(f"{name} must be a non-empty list of mappings.")
    if len(value) > _MAX_ENTRIES:
        raise ConfigurationError(f"{name} may contain at most {_MAX_ENTRIES} entries.")
    rows: list[dict[str, Any]] = []
    for index, item in enumerate(value):
        if not isinstance(item, Mapping):
            raise ConfigurationError(f"{name}[{index}] must be a mapping.")
        rows.append({str(key): item[key] for key in item})
    return rows


def _track_frame(value: object, name: str) -> Any:
    """Track metadata records as a DataFrame with validated ``name`` and ``strand`` columns."""

    pandas = require("pandas", _EXTRAS)
    rows = []
    for index, record in enumerate(_records_list(value, name)):
        track = _text(record.get("name"), f"{name}[{index}].name", maximum=256)
        strand = _choice(record.get("strand"), f"{name}[{index}].strand", _STRANDS)
        rows.append({**record, "name": track, "strand": strand})
    return pandas.DataFrame(rows)


def _anndata(value: object, name: str) -> Any:
    anndata = require("anndata", _EXTRAS)
    if not isinstance(value, anndata.AnnData):
        raise ConfigurationError(f"{name} must be an anndata.AnnData object.")
    return value


def _json_value(value: Any) -> Any:
    """JSON-safe value: NumPy and NaN become Python values; genome objects become text."""

    value = _clean(value)
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    return str(value)


def _records(frame: Any) -> list[dict[str, Any]]:
    pandas = require("pandas", _EXTRAS)
    return [
        {str(key): None if item is pandas.NA else _json_value(item) for key, item in row.items()}
        for row in frame.to_dict("records")
    ]


def _anndata_record(adata: Any) -> dict[str, Any]:
    return {
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "obs_names": [str(name) for name in adata.obs_names],
        "var_names": [str(name) for name in adata.var_names],
        "obs": _records(adata.obs),
        "var": _records(adata.var),
        "scores": _clean(adata.X),
        "layers": {str(name): _clean(adata.layers[name]) for name in sorted(adata.layers.keys())},
    }


def _variant_score_messages(value: object) -> list[Any]:
    atlas = require("alphagenome.protos.atlas_service_pb2", _EXTRAS)
    protobuf_message = require("google.protobuf.message", _EXTRAS)
    if isinstance(value, (str, bytes, bytearray)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError(
            "variant_scores must be a non-empty list of serialized DenseVariantScores bytes."
        )
    if len(value) > _MAX_ENTRIES:
        raise ConfigurationError(f"variant_scores may contain at most {_MAX_ENTRIES} entries.")
    messages = []
    for index, item in enumerate(value):
        if not isinstance(item, (bytes, bytearray)):
            raise ConfigurationError(f"variant_scores[{index}] must be bytes.")
        parsed = atlas.DenseVariantScores()
        try:
            parsed.ParseFromString(bytes(item))
        except protobuf_message.DecodeError as exc:
            raise _rejected(_ALPHAGENOME, exc) from exc
        messages.append(parsed)
    return messages


def _score_groups(value: object) -> list[list[Any]]:
    """Each AnnData is one group; a list of AnnData objects is one group."""

    anndata = require("anndata", _EXTRAS)
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence) or not value:
        raise ConfigurationError("scores must be a non-empty list of AnnData objects.")
    if len(value) > _MAX_ENTRIES:
        raise ConfigurationError(f"scores may contain at most {_MAX_ENTRIES} entries.")
    groups: list[list[Any]] = []
    for index, item in enumerate(value):
        if isinstance(item, anndata.AnnData):
            groups.append([item])
        elif (
            not isinstance(item, (str, bytes))
            and isinstance(item, Sequence)
            and item
            and all(isinstance(part, anndata.AnnData) for part in item)
        ):
            groups.append(list(item))
        else:
            raise ConfigurationError(
                f"scores[{index}] must be an AnnData object or a non-empty list of them."
            )
    return groups


def convert_alphagenome_variant_scores(
    variant_scores: Sequence[bytes],
    *,
    track_metadata: Mapping[str, Sequence[Mapping[str, Any]]],
    verbose: bool = True,
) -> dict[str, Any]:
    """Convert serialized ``DenseVariantScores`` messages to one AnnData per scorer.

    ``variant_scores`` holds ``DenseVariantScores.SerializeToString()`` bytes. ``track_metadata``
    maps each scorer name in the messages to its track records, because the protobuf payload
    does not always carry track names. The upstream converter is called without network access.
    """

    entries = _variant_score_messages(variant_scores)
    if not isinstance(track_metadata, Mapping) or not track_metadata:
        raise ConfigurationError("track_metadata must map each scorer name to track records.")
    pandas = require("pandas", _EXTRAS)
    frames: dict[str, Any] = {}
    for scorer, records in track_metadata.items():
        if not isinstance(scorer, str) or not scorer.strip():
            raise ConfigurationError("track_metadata keys must be non-empty scorer names.")
        frames[scorer.strip()] = pandas.DataFrame(
            _records_list(list(records), f"track_metadata[{scorer!r}]")
        )
    needed = {str(score.variant_scorer.name) for entry in entries for score in entry.scores}
    missing = sorted(needed - set(frames))
    if missing:
        raise ConfigurationError(f"track_metadata has no tracks for scorers: {missing}.")
    atlas = require("alphagenome.atlas.atlas", _EXTRAS)
    try:
        converted = atlas.convert_variant_scores_to_anndata(entries, frames)
    except (ValueError, KeyError, TypeError, IndexError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    return {
        "status": "completed",
        "operation": "convert_alphagenome_variant_scores",
        "variant_count": len(entries),
        "scorers": [
            {"scorer": name, "anndata": _anndata_record(converted[name])}
            for name in sorted(converted)
        ],
        "provenance": _provenance(
            _ALPHAGENOME, conversion="atlas.convert_variant_scores_to_anndata"
        ),
    }


def tidy_alphagenome_anndata(
    adata: Any,
    *,
    match_gene_strand: bool = True,
    include_extended_metadata: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Format one AlphaGenome score AnnData as a long table with one score per row.

    Rows whose gene strand and track strand disagree are removed when ``match_gene_strand``
    is true. ``include_extended_metadata`` adds biosample, assay and ontology columns.
    """

    scores = _anndata(adata, "adata")
    strand_match = _flag(match_gene_strand, "match_gene_strand")
    extended = _flag(include_extended_metadata, "include_extended_metadata")
    variant_scorers = require("alphagenome.models.v1.variant_scorers", _EXTRAS)
    try:
        frame = variant_scorers.tidy_anndata(
            scores,
            match_gene_strand=strand_match,
            include_extended_metadata=extended,
        )
    except (ValueError, KeyError, TypeError, AttributeError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    return {
        "status": "completed",
        "operation": "tidy_alphagenome_anndata",
        "row_count": len(frame),
        "columns": [str(column) for column in frame.columns],
        "records": _records(frame),
        "provenance": _provenance(_ALPHAGENOME, formatter="variant_scorers.tidy_anndata"),
    }


def tidy_alphagenome_scores(
    scores: Sequence[Any],
    *,
    match_gene_strand: bool = True,
    include_extended_metadata: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Format several AlphaGenome score AnnData objects as one long table.

    ``scores`` is a list of AnnData objects, or a list of lists of AnnData objects (the shape
    returned for several variants and scorers). Results keep the input order.
    """

    groups = _score_groups(scores)
    strand_match = _flag(match_gene_strand, "match_gene_strand")
    extended = _flag(include_extended_metadata, "include_extended_metadata")
    variant_scorers = require("alphagenome.models.v1.variant_scorers", _EXTRAS)
    try:
        frame = variant_scorers.tidy_scores(
            groups,
            match_gene_strand=strand_match,
            include_extended_metadata=extended,
        )
    except (ValueError, KeyError, TypeError, AttributeError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    if frame is None:
        frame_records: list[dict[str, Any]] = []
        columns: list[str] = []
    else:
        frame_records = _records(frame)
        columns = [str(column) for column in frame.columns]
    return {
        "status": "completed",
        "operation": "tidy_alphagenome_scores",
        "group_count": len(groups),
        "row_count": len(frame_records),
        "columns": columns,
        "records": frame_records,
        "provenance": _provenance(_ALPHAGENOME, formatter="variant_scorers.tidy_scores"),
    }


def merge_alphagenome_track_metadata(
    track_metadata: Sequence[Mapping[str, Any]],
    *,
    validate_tracks: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge stranded track records (``+`` and ``-``) into single unstranded tracks.

    Each record needs ``name`` and ``strand`` (``+``, ``-`` or ``.``). Unstranded records with
    the same name as a merged track are dropped. With ``validate_tracks``, the positive and
    negative tracks must have the same ordered names.
    """

    frame = _track_frame(track_metadata, "track_metadata")
    validate = _flag(validate_tracks, "validate_tracks")
    utils = require("alphagenome.models.v1.variant_scoring_utils", _EXTRAS)
    try:
        merged = utils.merge_stranded_track_metadata(frame, validate_tracks=validate)
    except (ValueError, KeyError, TypeError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    return {
        "status": "completed",
        "operation": "merge_alphagenome_track_metadata",
        "track_count": len(merged),
        "tracks": _records(merged),
        "provenance": _provenance(_ALPHAGENOME, merge="variant_scoring_utils"),
    }


def merge_alphagenome_stranded_gene_tracks(
    scores: Any,
    *,
    validate_tracks: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Merge gene scores from stranded tracks into one unstranded track per name.

    For each gene, the positive track is used when the gene is on ``+`` and the negative
    track when it is on ``-``. If the AnnData has no ``gene_id`` and ``strand`` columns in
    ``obs``, the scores are returned unchanged.
    """

    adata = _anndata(scores, "scores")
    validate = _flag(validate_tracks, "validate_tracks")
    utils = require("alphagenome.models.v1.variant_scoring_utils", _EXTRAS)
    try:
        merged = utils.merge_stranded_gene_tracks(adata, validate_tracks=validate)
    except (ValueError, KeyError, TypeError, IndexError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    return {
        "status": "completed",
        "operation": "merge_alphagenome_stranded_gene_tracks",
        "anndata": _anndata_record(merged),
        "provenance": _provenance(_ALPHAGENOME, merge="variant_scoring_utils"),
    }


def unmerge_alphagenome_stranded_gene_tracks(
    scores: Any,
    *,
    track_metadata: Sequence[Mapping[str, Any]],
    validate_tracks: bool = True,
    verbose: bool = True,
) -> dict[str, Any]:
    """Split unstranded gene scores back into stranded tracks using the original metadata.

    Each gene keeps the track matching its strand; the mismatched strand columns are NaN
    (``null`` in the output). ``track_metadata`` needs ``name`` and ``strand`` per track.
    """

    adata = _anndata(scores, "scores")
    frame = _track_frame(track_metadata, "track_metadata")
    validate = _flag(validate_tracks, "validate_tracks")
    utils = require("alphagenome.models.v1.variant_scoring_utils", _EXTRAS)
    try:
        unmerged = utils.unmerge_stranded_gene_tracks(
            adata, track_metadata=frame, validate_tracks=validate
        )
    except (ValueError, KeyError, TypeError, IndexError) as exc:
        raise _rejected(_ALPHAGENOME, exc) from exc
    return {
        "status": "completed",
        "operation": "unmerge_alphagenome_stranded_gene_tracks",
        "anndata": _anndata_record(unmerged),
        "provenance": _provenance(_ALPHAGENOME, unmerge="variant_scoring_utils"),
    }


def list_alphagenome_recommended_scorers(
    organism: str = "homo_sapiens",
    *,
    verbose: bool = True,
) -> dict[str, Any]:
    """List the recommended variant scorers for an organism, in name order.

    ``organism`` is ``homo_sapiens`` or ``mus_musculus``. Each entry gives the registry name
    and the scorer description.
    """

    name = _choice(organism, "organism", _ORGANISMS)
    model = require("alphagenome.models.v1.dna_model", _EXTRAS)
    variant_scorers = require("alphagenome.models.v1.variant_scorers", _EXTRAS)
    species = model.Organism[name.upper()].value
    chosen = {id(scorer) for scorer in variant_scorers.get_recommended_scorers(species)}
    registry = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    scorers = [
        {"name": key, "scorer": str(registry[key])}
        for key in sorted(registry)
        if id(registry[key]) in chosen
    ]
    return {
        "status": "completed",
        "operation": "list_alphagenome_recommended_scorers",
        "organism": name,
        "scorer_count": len(scorers),
        "scorers": scorers,
        "provenance": _provenance(
            _ALPHAGENOME, selection="variant_scorers.get_recommended_scorers"
        ),
    }


__all__ = [
    "convert_alphagenome_variant_scores",
    "list_alphagenome_recommended_scorers",
    "merge_alphagenome_stranded_gene_tracks",
    "merge_alphagenome_track_metadata",
    "tidy_alphagenome_anndata",
    "tidy_alphagenome_scores",
    "unmerge_alphagenome_stranded_gene_tracks",
]

_register_public_verbose(globals(), tuple(__all__))
